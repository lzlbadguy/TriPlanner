// mvie2d_pentagon.hpp
// -----------------------------------------------------------------------------
// Header-only C++17 library for Inscribed Ellipse in a Convex Pentagon (N=5)
// Method: Linear constraints in the dual line-space (one-shot linear algebra).
// Inputs supported:
//   - 5 Halfspaces (a^T x <= b)
//   - 5 Tangency lines (ax + by + c <= 0)
//   - 5 CCW vertices (convex pentagon)
// Output: mvie2d::Ellipse2  (same semantics as menn_tri_se / menn_quad_se):
//   (x - c_center)^T X (x - c_center) <= 1
// Requires: Eigen >= 3.3
// License: MIT
// -----------------------------------------------------------------------------

#pragma once
#include <Eigen/Dense>
#include <array>
#include <vector>
#include <optional>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

#include "mvie2d_types.hpp" // ← 统一 TangencyLine / Ellipse2 / Line2

namespace mvie2d
{

    // ----------------- 工具结构与函数 -----------------

    inline Eigen::Vector2d rotate90R(const Eigen::Vector2d &v)
    {
        return {v.y(), -v.x()};
    }

    // 顶点(CCW) → 外法向 halfspaces（a^T x <= b）
    inline std::vector<Line2> buildHalfspacesFromVerticesCCW(const std::vector<Eigen::Vector2d> &V)
    {
        const int n = (int)V.size();
        std::vector<Line2> H;
        H.reserve(n);
        for (int i = 0; i < n; ++i)
        {
            Eigen::Vector2d p = V[i], q = V[(i + 1) % n];
            Eigen::Vector2d e = q - p;
            if (e.norm() < 1e-15)
                continue;
            // CCW 顶点序：右旋得外法向
            Eigen::Vector2d n_out = rotate90R(e);
            const double len = n_out.norm();
            if (len < 1e-15)
                continue;
            n_out /= len;
            double b = n_out.dot(p);
            H.push_back({n_out, b});
        }
        return H;
    }

    // ax+by+c <= 0 → halfspace n·x <= b，n=(a,b)/|| (a,b) ||，b = -c/|| (a,b) ||
    template <class TangencyLike>
    inline std::vector<Line2> halfspacesFromABC5(const std::vector<TangencyLike> &TL)
    {
        std::vector<Line2> H;
        H.reserve(TL.size());
        for (const auto &t : TL)
        {
            Eigen::Vector2d n(t.a, t.b);
            const double norm = n.norm();
            if (norm < 1e-15)
                continue; // 跳过退化
            n /= norm;
            const double bb = -t.c / norm;
            H.push_back({n, bb});
        }
        return H;
    }

    // halfspace -> 齐次直线（线域）：l = [a_x, a_y, -b]
    inline Eigen::Vector3d line_to_homog(const Line2 &h)
    {
        return Eigen::Vector3d(h.a.x(), h.a.y(), -h.b);
    }

    // 根据 X 填充椭圆的派生量（a,b,theta）
    inline void fill_derived_from_X(Ellipse2 &E)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(E.X);
        if (es.info() != Eigen::Success)
        {
            E.valid = false;
            E.major_axis = E.minor_axis = 0.0;
            E.c_theta = 0.0;
            return;
        }

        auto lam = es.eigenvalues(); // 升序
        auto V = es.eigenvectors();

        if (lam(0) <= 0.0 || lam(1) <= 0.0)
        {
            E.valid = false;
            E.major_axis = E.minor_axis = 0.0;
            E.c_theta = 0.0;
            return;
        }

        double a = 1.0 / std::sqrt(lam(0));
        double b = 1.0 / std::sqrt(lam(1));
        Eigen::Vector2d v = V.col(0);
        double theta = std::atan2(v.y(), v.x());

        // 保证 major_axis >= minor_axis
        if (b > a)
        {
            std::swap(a, b);
            v = V.col(1);
            theta = std::atan2(v.y(), v.x());
        }

        E.major_axis = a;
        E.minor_axis = b;
        E.c_theta = theta;
    }

    // 点域二次型 MP => 椭圆（含导出量）
    // MP = [A b; b^T c], 满足 x^T A x + 2 b^T x + c = 0
    // 推导出：(x - center)^T Q (x - center) = 1, 其中 Q = A / k', k' = b^T A^{-1} b - c
    inline std::optional<Ellipse2> ellipse_from_Mp(const Eigen::Matrix3d &MP)
    {
        Eigen::Matrix2d A;
        A << MP(0, 0), MP(0, 1),
            MP(1, 0), MP(1, 1);
        Eigen::Vector2d b(MP(0, 2), MP(1, 2));
        double c = MP(2, 2);

        if (std::abs(A.determinant()) < 1e-18)
            return std::nullopt;

        Eigen::Matrix2d Ainv = A.inverse();
        Eigen::Vector2d center = -Ainv * b;
        double kprime = b.transpose() * Ainv * b - c; // > 0 for ellipse
        if (!(kprime > 0.0))
            return std::nullopt;

        Eigen::Matrix2d Q = A / kprime; // 应为 SPD
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> esQ(Q);
        if (esQ.eigenvalues().minCoeff() <= 0.0)
            return std::nullopt;

        Ellipse2 E;
        E.c_center = center;
        E.X = Q; // 统一语义：(x-c)^T X (x-c) <= 1
        E.valid = true;

        fill_derived_from_X(E); // 填充 a,b,theta
        if (!E.valid)
            return std::nullopt;

        return E;
    }

    // ----------------- 五边形求解（核心） -----------------
    // 输入：5 条 halfspaces（a^T x <= b）
    // 方法：在“线域”构造 5x6 线性约束，解其零空间（1维，唯一到尺度）得到对偶二次型 ML，
    //       再取点域 MP = ML^{-1} 并恢复椭圆。
    inline std::optional<Ellipse2> inscribedEllipsePentagon(const std::array<Line2, 5> &H5)
    {
        Eigen::Matrix<double, 5, 6> A; // 行 i 对应第 i 条边的约束
        for (int i = 0; i < 5; ++i)
        {
            Eigen::Vector3d l = line_to_homog(H5[i]); // l = [ax, ay, c]
            const double ax = l(0), ay = l(1), c = l(2);
            A(i, 0) = ax * ax;     // u1
            A(i, 1) = 2 * ax * ay; // u2
            A(i, 2) = 2 * ax * c;  // u3
            A(i, 3) = ay * ay;     // u4
            A(i, 4) = 2 * ay * c;  // u5
            A(i, 5) = c * c;       // u6
        }

        // SVD 求零空间（到尺度唯一） → ML 对偶二次型
        Eigen::JacobiSVD<Eigen::Matrix<double, 5, 6>> svd(A, Eigen::ComputeFullV);
        Eigen::Matrix<double, 6, 1> m = svd.matrixV().col(5); // 最小奇值对应列

        Eigen::Matrix3d ML;
        ML << m(0), m(1), m(2),
            m(1), m(3), m(4),
            m(2), m(4), m(5);

        if (std::abs(ML.determinant()) < 1e-18)
            return std::nullopt;

        // 点域二次型
        Eigen::Matrix3d MP = ML.inverse();

        // 先试正常符号
        if (auto E = ellipse_from_Mp(MP))
        {
            return E;
        }
        // 若失败，翻转整体符号再试一次（零空间向量到尺度唯一，±等价）
        return ellipse_from_Mp(-MP);
    }

    inline Ellipse2 inscribedEllipsePentagon(const std::vector<TangencyLine> &TL)
    {
        Ellipse2 bad;
        if (TL.size() != 5)
            return bad;

        auto H = halfspacesFromABC5(TL);
        if (H.size() != 5)
            return bad;

        std::array<Line2, 5> H5;
        for (int i = 0; i < 5; i++)
            H5[i] = H[i];

        auto Eopt = inscribedEllipsePentagon(H5); // 内部还是 optional
        if (!Eopt)
            return bad;
        return *Eopt;
    }

    // 重载 2：5 个 CCW 顶点（要求凸）
    inline std::optional<Ellipse2> inscribedEllipsePentagon(const std::vector<Eigen::Vector2d> &verticesCCW)
    {
        if (verticesCCW.size() != 5)
            return std::nullopt;
        auto H = buildHalfspacesFromVerticesCCW(verticesCCW);
        if (H.size() != 5)
            return std::nullopt;
        std::array<Line2, 5> H5;
        for (int i = 0; i < 5; ++i)
            H5[i] = H[i];
        return inscribedEllipsePentagon(H5);
    }

} // namespace mvie2d
