// ============================================================================
// menn_tri_se.hpp  — Steiner inellipse (MENN of a triangle)
// From 3 tangency lines: a x + b y + c <= 0  (each line is a polygon edge)
// Output style mirrors menn_quad_se.hpp / menn_pentagon.hpp:
//   struct Ellipse2 { c_center, X, valid, major_axis, minor_axis, c_theta, area() }
// License: MIT
// ============================================================================

#pragma once
#include <Eigen/Dense>
#include <array>
#include <vector>
#include <optional>
#include <cmath>
#include <limits>
#include <algorithm>

#include "mvie2d_types.hpp" // ← 统一的 TangencyLine / Ellipse2 定义

namespace mvie2d
{

    namespace detail
    {

        inline std::optional<Eigen::Vector2d>
        intersect(double a1, double b1, double c1,
                  double a2, double b2, double c2,
                  double eps = 1e-12)
        {
            // a1 x + b1 y + c1 = 0
            // a2 x + b2 y + c2 = 0
            const double det = a1 * b2 - a2 * b1;
            if (std::abs(det) < eps)
                return std::nullopt;
            Eigen::Vector2d p;
            p.x() = (b1 * c2 - b2 * c1) / det;
            p.y() = (a2 * c1 - a1 * c2) / det;
            return p;
        }

        inline void fill_derived(Ellipse2 &E)
        {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(E.X);
            if (es.info() != Eigen::Success)
            {
                E.valid = false;
                return;
            }

            auto lam = es.eigenvalues(); // ascending
            auto V = es.eigenvectors();

            if (lam(0) <= 0.0 || lam(1) <= 0.0)
            {
                E.valid = false;
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

    } // namespace detail

    // ---- main solver: from 3 tangency lines (triangle) ----
    // TL.size() must be 3; inequalities must all face "inside": a x + b y + c <= 0
    inline std::optional<Ellipse2>
    inscribedEllipseTriangle_ori(const std::vector<TangencyLine> &TL, double eps = 1e-12)
    {
        if (TL.size() != 3)
            return std::nullopt;

        // Vertex convention (keeps Steiner formulas consistent):
        //   v1 = L2 ∩ L3,  v2 = L3 ∩ L1,  v3 = L1 ∩ L2
        auto v1 = detail::intersect(TL[1].a, TL[1].b, TL[1].c,
                                    TL[2].a, TL[2].b, TL[2].c, eps);
        auto v2 = detail::intersect(TL[2].a, TL[2].b, TL[2].c,
                                    TL[0].a, TL[0].b, TL[0].c, eps);
        auto v3 = detail::intersect(TL[0].a, TL[0].b, TL[0].c,
                                    TL[1].a, TL[1].b, TL[1].c, eps);

        Ellipse2 E;
        if (!v1 || !v2 || !v3)
        {
            E.valid = false;
            return E;
        }

        const Eigen::Vector2d V1 = *v1;
        const Eigen::Vector2d V2 = *v2;
        const Eigen::Vector2d V3 = *v3;

        // Steiner inellipse for triangle:
        // center v0 = (V1+V2+V3)/3
        // f1 = 1/2 (v0 - V3)
        // f2 = 1/(2√3) (V1 - V2)
        const Eigen::Vector2d v0 = (V1 + V2 + V3) / 3.0;
        const Eigen::Vector2d f1 = 0.5 * (v0 - V3);
        const Eigen::Vector2d f2 = (1.0 / (2.0 * std::sqrt(3.0))) * (V1 - V2);

        Eigen::Matrix2d F;
        F.col(0) = f1;
        F.col(1) = f2;
        Eigen::Matrix2d S = F * F.transpose(); // SPD for non-degenerate triangle
        if (S.determinant() <= eps)
        {
            E.valid = false;
            return E;
        }

        E.c_center = v0;
        E.X = S.inverse();
        E.valid = true;

        detail::fill_derived(E);
        return E;
    }

    inline Ellipse2
    inscribedEllipseTriangle(const std::vector<TangencyLine> &TL, double eps = 1e-12)
    {
        Ellipse2 bad;
        bad.valid = false;

        auto Eopt = inscribedEllipseTriangle_ori(TL, eps);
        if (!Eopt)
            return bad;
        return *Eopt;
    }

} // namespace mvie2d
