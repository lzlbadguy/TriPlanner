#pragma once
#include <Eigen/Dense>
#include <cmath>

namespace mvie2d
{

    // ----------------- Halfspace & Tangency -----------------

    // halfspace: a^T x <= b
    struct Line2
    {
        Eigen::Vector2d a{0.0, 0.0}; // 法向量（不一定单位）
        double b{0.0};

        Line2() = default;
        Line2(const Eigen::Vector2d &n, double bb) : a(n), b(bb) {}
        Line2(double ax, double ay, double bb) : a(ax, ay), b(bb) {}
    };

    // ax + by + c <= 0 表示“内侧”
    struct TangencyLine
    {
        double a{0.0}, b{0.0}, c{0.0};
        TangencyLine() = default;
        TangencyLine(double A, double B, double C) : a(A), b(B), c(C) {}
    };

    // ----------------- Ellipse -----------------

    // 椭圆用 (x-c)^T X (x-c) <= 1 统一表示（X 为 SPD）
    struct Ellipse2
    {
        bool valid{false};
        Eigen::Vector2d c_center{0.0, 0.0};
        Eigen::Matrix2d X = Eigen::Matrix2d::Identity(); // SPD

        // 派生量（可选，某些库会填充它们）
        double major_axis{0.0}; // a
        double minor_axis{0.0}; // b
        double c_theta{0.0};    // 主轴角

        inline double area() const
        {
            if (!valid)
                return 0.0;
            const double detX = X.determinant();
            if (!(detX > 0.0))
                return 0.0;
            // 对于 (x-c)^T X (x-c) <= 1, area = π / sqrt(det(X)) = π a b
            return M_PI / std::sqrt(detX);
        }
    };

    // —— 兼容辅助（给旧代码用：需要 L 或 c 时怎么取） ——
    inline Eigen::Matrix2d choleskyL(const Ellipse2 &E)
    {
        Eigen::LLT<Eigen::Matrix2d> llt(E.X);
        return llt.matrixL(); // X = L L^T
    }

    inline Eigen::Vector2d &center(Ellipse2 &E) { return E.c_center; }
    inline const Eigen::Vector2d &center(const Ellipse2 &E) { return E.c_center; }

} // namespace mvie2d
