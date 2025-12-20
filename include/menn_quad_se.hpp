#pragma once
#include <Eigen/Dense>
#include <array>
#include <vector>
#include <optional>
#include <limits>
#include <cmath>
#include <algorithm>

#include "mvie2d_types.hpp" // ← 统一的 TangencyLine / Ellipse2 定义

// =====================================================================================
//  Ellipse of MAXIMAL AREA inscribed in a convex quadrilateral (Horwitz 2005)
//  Core solver works in a normalized coord. system, then maps back.
//  输出：mvie2d::Ellipse2，语义与 menn_tri_se.hpp 一致：
//    (x - c_center)^T X (x - c_center) <= 1
// =====================================================================================

namespace mvie2d
{
    struct Options
    {
        double eps = 1e-12;
    };

    // --------------------------------- Utils ------------------------------------

    inline Eigen::Vector2d rot90L(const Eigen::Vector2d &v) { return {-v.y(), v.x()}; }

    inline bool line_intersection(const Line2 &l1, const Line2 &l2, Eigen::Vector2d &x, double eps = 1e-12)
    {
        double det = l1.a.x() * l2.a.y() - l1.a.y() * l2.a.x();
        if (std::abs(det) < eps)
            return false;
        x.x() = (l1.b * l2.a.y() - l1.a.y() * l2.b) / det;
        x.y() = (-l1.b * l2.a.x() + l1.a.x() * l2.b) / det;
        return x.allFinite();
    }

    inline bool inside_all(const std::vector<Line2> &H, const Eigen::Vector2d &p, double eps)
    {
        for (const auto &h : H)
            if (h.a.dot(p) - h.b > eps)
                return false;
        return true;
    }

    inline std::vector<Eigen::Vector2d>
    quad_vertices_from_halfspaces(const std::vector<Line2> &H, double eps)
    {
        std::vector<Eigen::Vector2d> pts;
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
            {
                Eigen::Vector2d p;
                if (line_intersection(H[i], H[j], p, eps) && inside_all(H, p, 1e-9))
                    pts.push_back(p);
            }
        // unique filter
        std::vector<Eigen::Vector2d> uniq;
        const double tol2 = 1e-16;
        for (auto &p : pts)
        {
            bool ok = true;
            for (auto &q : uniq)
                if ((p - q).squaredNorm() < tol2)
                {
                    ok = false;
                    break;
                }
            if (ok)
                uniq.push_back(p);
        }
        // centroid
        Eigen::Vector2d c(0, 0);
        for (auto &p : uniq)
            c += p;
        if (!uniq.empty())
            c /= double(uniq.size());
        // angle sort
        std::sort(uniq.begin(), uniq.end(), [&](const Eigen::Vector2d &p, const Eigen::Vector2d &q)
                  { return std::atan2(p.y() - c.y(), p.x() - c.x()) < std::atan2(q.y() - c.y(), q.x() - c.x()); });
        return uniq;
    }

    struct Affine2
    {
        Eigen::Matrix2d A;
        Eigen::Vector2d t;
    };

    inline Affine2 make_normalization(const std::array<Eigen::Vector2d, 4> &V)
    {
        Eigen::Matrix2d B;
        B.col(0) = V[1] - V[0];
        B.col(1) = V[2] - V[0];
        Eigen::Matrix2d A = B.inverse();
        Eigen::Vector2d t = -A * V[0];
        return {A, t};
    }

    inline Eigen::Vector2d apply_affine(const Affine2 &T, const Eigen::Vector2d &x) { return T.A * x + T.t; }
    inline Eigen::Vector2d apply_affine_inv(const Affine2 &T, const Eigen::Vector2d &x) { return T.A.inverse() * (x - T.t); }

    // ------------------------------ Core (normalized) --------------------------

    struct EllipseNormalized
    {
        Eigen::Vector2d c; // center in normalized coords
        Eigen::Matrix2d X; // 在归一化坐标系下的“轴长平方矩阵” M，特征值 ~ a^2,b^2
        bool valid{false};
    };

    // 给定中心 c 以及三条不共线的边（ax + by + c = 0），解出对称矩阵 X（此处 X=M）
    static inline std::optional<Eigen::Matrix2d>
    solve_X_from_three_lines(const Eigen::Vector2d &c,
                             const std::array<Eigen::Vector3d, 3> &abc)
    {
        auto unit_from_abc = [](double a, double b, double cc,
                                Eigen::Vector2d &n, double &d)
        {
            double len = std::hypot(a, b);
            if (len < 1e-15)
            {
                n = Eigen::Vector2d::Zero();
                d = 0.0;
                return;
            }
            n = Eigen::Vector2d(a / len, b / len);
            d = -cc / len;
        };
        auto row_from_n = [](const Eigen::Vector2d &n)
        { return Eigen::RowVector3d(n.x() * n.x(), 2.0 * n.x() * n.y(), n.y() * n.y()); };
        auto rhs_from_nd = [&](const Eigen::Vector2d &n, double d)
        { double v = n.dot(c) - d; return v*v; };

        Eigen::Vector2d n1, n2, n3;
        double d1 = 0, d2 = 0, d3 = 0;
        // 逐条线取 a,b,c 分量
        unit_from_abc(abc[0].x(), abc[0].y(), abc[0].z(), n1, d1);
        unit_from_abc(abc[1].x(), abc[1].y(), abc[1].z(), n2, d2);
        unit_from_abc(abc[2].x(), abc[2].y(), abc[2].z(), n3, d3);

        Eigen::Matrix3d A3;
        A3.row(0) = row_from_n(n1);
        A3.row(1) = row_from_n(n2);
        A3.row(2) = row_from_n(n3);
        Eigen::Vector3d b3;
        b3 << rhs_from_nd(n1, d1), rhs_from_nd(n2, d2), rhs_from_nd(n3, d3);

        if (A3.fullPivLu().rank() < 3)
            return std::nullopt;

        Eigen::Vector3d xvec = A3.fullPivLu().solve(b3);
        Eigen::Matrix2d X;
        X << xvec(0), xvec(1), xvec(1), xvec(2);
        return X; // 这里的 X 就是 M（特征值 a^2,b^2）
    }

    std::optional<EllipseNormalized>
    horwitz_max_area_normalized(double s, double t, const Options &opt = {})
    {
        // 1. solve for h*,k*
        auto L_of_h = [&](double h)
        {
            return ((s - t) + 2.0 * h * (t - 1.0)) / (2.0 * (s - 1.0));
        };
        double h_lo = 0.5, h_hi = 0.5 * s;

        if (h_hi < h_lo)
        {
            std::swap(h_lo, h_hi);
        }

        auto clampI = [&](double h)
        { return std::clamp(h, h_lo + 1e-12, h_hi - 1e-12); };

        double A = 12.0 * (1.0 - t);
        double B = 4.0 * (s * t - 2.0 * s + t - 1.0);
        double C = s * s - s * t + 2.0 * s;
        double h_star;

        if (std::abs(A) < 1e-14)
        {
            if (std::abs(B) < 1e-18)
                return std::nullopt;
            h_star = clampI(-C / B);
        }
        else
        {
            double D = B * B - 4 * A * C;
            if (D < 0)
            {
                return std::nullopt;
            }
            double sqrtD = std::sqrt(std::max(0.0, D));
            double q = -0.5 * (B + std::copysign(sqrtD, B));
            double r1 = q / A, r2 = (std::abs(q) < 1e-18) ? r1 : C / q;
            bool r1_in = (r1 > h_lo && r1 < h_hi), r2_in = (r2 > h_lo && r2 < h_hi);
            h_star = r1_in && !r2_in ? r1 : (r2_in && !r1_in ? r2 : (std::abs(r1 - 0.5 * (h_lo + h_hi)) < std::abs(r2 - 0.5 * (h_lo + h_hi)) ? r1 : r2));
            h_star = clampI(h_star);
        }
        double k_star = L_of_h(h_star);
        Eigen::Vector2d c(h_star, k_star);

        // 2. build 3x3 system (四条边中的三条)
        auto unit_from_abc = [](double a, double b, double cc, Eigen::Vector2d &n, double &d)
        {
            double len = std::hypot(a, b);
            n = Eigen::Vector2d(a / len, b / len);
            d = -cc / len;
        };

        Eigen::Vector2d n1, n2, n3, n4;
        double d1, d2, d3, d4;
        unit_from_abc(0, 1, 0, n1, d1);
        unit_from_abc(1, 0, 0, n2, d2);
        unit_from_abc(t, -(s - 1.0), -t, n3, d3);
        unit_from_abc((t - 1.0), -s, s, n4, d4);

        auto row_from_n = [](const Eigen::Vector2d &n)
        { return Eigen::RowVector3d(n.x() * n.x(), 2 * n.x() * n.y(), n.y() * n.y()); };
        auto rhs_from_nd = [&](const Eigen::Vector2d &n, double d)
        {double val=n.dot(c)-d;return val*val; };

        Eigen::Matrix3d A3;
        A3.row(0) = row_from_n(n1);
        A3.row(1) = row_from_n(n2);
        A3.row(2) = row_from_n(n3);
        Eigen::Vector3d b3;
        b3 << rhs_from_nd(n1, d1), rhs_from_nd(n2, d2), rhs_from_nd(n3, d3);

        if (A3.fullPivLu().rank() < 3)
            return std::nullopt;
        Eigen::Vector3d xvec = A3.fullPivLu().solve(b3);
        Eigen::Matrix2d X;
        X << xvec(0), xvec(1), xvec(1), xvec(2);

        return EllipseNormalized{c, X, true}; // 这里的 X 仍然是 M
    }

    // ------------------------------ Public API ---------------------------------

    inline std::optional<Ellipse2>
    computeMaxAreaInscribedEllipse(const std::vector<Line2> &H, const Options &opt = {})
    {
        if (H.size() != 4)
            return std::nullopt;
        auto Vv = quad_vertices_from_halfspaces(H, opt.eps);
        if (Vv.size() != 4)
            return std::nullopt;
        std::array<Eigen::Vector2d, 4> V{Vv[0], Vv[1], Vv[3], Vv[2]};

        auto T = make_normalization(V);
        Eigen::Vector2d v3n = apply_affine(T, V[3]);
        double s = v3n.x(), t = v3n.y();

        auto is_parallelogram_st = [&](double s, double t) -> bool
        {
            const double tol = 1e-10;
            return ((Eigen::Vector2d(s, t) - Eigen::Vector2d(1.0, 1.0)).norm() <= tol);
        };

        // --------- Case 1：平行四边形（Steiner 椭圆）---------
        if (is_parallelogram_st(s, t))
        {
            Eigen::Vector2d c = (V[0] + V[1] + V[2] + V[3]) / 4.0;
            Eigen::Vector2d e1 = V[1] - V[0];
            Eigen::Vector2d e2 = V[2] - V[0];
            Eigen::Matrix2d B;
            B.col(0) = e1;
            B.col(1) = e2;

            // 这里的 M 是 “轴长平方矩阵”
            Eigen::Matrix2d M = 0.25 * (B * B.transpose());

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(M);
            if (es.info() != Eigen::Success)
                return std::nullopt;
            auto eval = es.eigenvalues();
            auto evec = es.eigenvectors();
            if (eval.minCoeff() <= 0)
                return std::nullopt;

            int i_max = (eval(1) >= eval(0)) ? 1 : 0;
            int i_min = 1 - i_max;
            double a = std::sqrt(eval(i_max));
            double b = std::sqrt(eval(i_min));
            Eigen::Vector2d u_major = evec.col(i_max);

            Ellipse2 E;
            E.valid = true;
            E.c_center = c;
            E.major_axis = a;
            E.minor_axis = b;
            E.c_theta = std::atan2(u_major.y(), u_major.x());
            E.X = M.inverse(); // 统一语义：(x-c)^T X (x-c) <= 1

            return E;
        }

        const double tol = 1e-10;

        // --------- Case 2：梯形 A（竖边平行，s≈1）---------
        if (std::abs(s - 1.0) <= tol)
        {
            Eigen::Vector2d c_n(0.5, (t + 1.0) * 0.25);

            std::array<Eigen::Vector3d, 3> abc = {
                Eigen::Vector3d(0.0, 1.0, 0.0),        // y = 0
                Eigen::Vector3d(1.0, 0.0, -1.0),       // x = 1
                Eigen::Vector3d(-(t - 1.0), 1.0, -1.0) // y - 1 - (t-1)x = 0
            };

            auto Xn_opt = solve_X_from_three_lines(c_n, abc);
            if (!Xn_opt)
                return std::nullopt;

            Eigen::Matrix2d Ainv = T.A.inverse();
            Eigen::Vector2d c_orig = apply_affine_inv(T, c_n);
            Eigen::Matrix2d M_orig = Ainv * (*Xn_opt) * Ainv.transpose(); // M

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(M_orig);
            if (es.info() != Eigen::Success || es.eigenvalues().minCoeff() <= 0)
                return std::nullopt;

            auto eval = es.eigenvalues();
            auto evec = es.eigenvectors();
            int i_max = (eval(1) >= eval(0)) ? 1 : 0;
            int i_min = 1 - i_max;
            double a = std::sqrt(eval(i_max));
            double b = std::sqrt(eval(i_min));
            Eigen::Vector2d u_major = evec.col(i_max);

            Ellipse2 E;
            E.valid = true;
            E.c_center = c_orig;
            E.major_axis = a;
            E.minor_axis = b;
            E.c_theta = std::atan2(u_major.y(), u_major.x());
            E.X = M_orig.inverse();

            return E;
        }

        // --------- Case 3：梯形 B（横边平行，t≈1）---------
        if (std::abs(t - 1.0) <= tol)
        {
            Eigen::Vector2d c_n((s + 1.0) * 0.25, 0.5);

            std::array<Eigen::Vector3d, 3> abc = {
                Eigen::Vector3d(1.0, 0.0, 0.0),       // x = 0
                Eigen::Vector3d(0.0, 1.0, 0.0),       // y = 0
                Eigen::Vector3d(-1.0, (s - 1.0), 1.0) // (s-1)y - x + 1 = 0
            };

            auto Xn_opt = solve_X_from_three_lines(c_n, abc);
            if (!Xn_opt)
                return std::nullopt;

            Eigen::Matrix2d Ainv = T.A.inverse();
            Eigen::Vector2d c_orig = apply_affine_inv(T, c_n);
            Eigen::Matrix2d M_orig = Ainv * (*Xn_opt) * Ainv.transpose();

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(M_orig);
            if (es.info() != Eigen::Success || es.eigenvalues().minCoeff() <= 0)
                return std::nullopt;

            auto eval = es.eigenvalues();
            auto evec = es.eigenvectors();
            int i_max = (eval(1) >= eval(0)) ? 1 : 0;
            int i_min = 1 - i_max;
            double a = std::sqrt(eval(i_max));
            double b = std::sqrt(eval(i_min));
            Eigen::Vector2d u_major = evec.col(i_max);

            Ellipse2 E;
            E.valid = true;
            E.c_center = c_orig;
            E.major_axis = a;
            E.minor_axis = b;
            E.c_theta = std::atan2(u_major.y(), u_major.x());
            E.X = M_orig.inverse();

            return E;
        }

        // --------- Case 4：一般四边形 ----------
        auto En = horwitz_max_area_normalized(s, t, opt);
        if (!En || !En->valid)
            return std::nullopt;

        Eigen::Matrix2d Ainv = T.A.inverse();
        Eigen::Vector2d c_orig = apply_affine_inv(T, En->c);
        Eigen::Matrix2d M_orig = Ainv * En->X * Ainv.transpose(); // M

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(M_orig);
        if (es.info() != Eigen::Success)
            return std::nullopt;
        auto eval = es.eigenvalues();
        auto evec = es.eigenvectors();
        if (eval.minCoeff() <= 0)
            return std::nullopt;

        int i_max = (eval(1) >= eval(0)) ? 1 : 0;
        int i_min = 1 - i_max;
        double a = std::sqrt(eval(i_max));
        double b = std::sqrt(eval(i_min));
        Eigen::Vector2d u_major = evec.col(i_max);

        Ellipse2 E;
        E.valid = true;
        E.c_center = c_orig;
        E.major_axis = a;
        E.minor_axis = b;
        E.c_theta = std::atan2(u_major.y(), u_major.x());
        E.X = M_orig.inverse(); // 输出统一为 (x-c)^T X (x-c) <= 1

        return E;
    }

    // -------------------------- From Tangency Lines ----------------------------

    template <class TangencyLike>
    inline std::vector<Line2>
    halfspacesFromABC4(const std::vector<TangencyLike> &TL)
    {
        std::vector<Line2> H;
        H.reserve(TL.size());
        for (const auto &t : TL)
        {
            Eigen::Vector2d n(t.a, t.b);
            double c = t.c;
            double nn = n.norm();
            if (nn < 1e-15)
                continue;
            n /= nn;
            double bb = -c / nn;
            H.push_back({n, bb});
        }
        return H;
    }

    template <class TangencyLike>
    inline Ellipse2 inscribedEllipseQuad(const std::vector<TangencyLike> &TL, const Options &opt = {})
    {
        auto H = halfspacesFromABC4(TL);
        auto Eopt = computeMaxAreaInscribedEllipse(H, opt);
        if (Eopt)
            return *Eopt;
        Ellipse2 bad;
        bad.valid = false;
        return bad;
    }
} // namespace mvie2d
