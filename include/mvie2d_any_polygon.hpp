#pragma once
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <optional>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "mvie2d_types.hpp"  // Line2, TangencyLine, Ellipse2
#include "menn_tri_se.hpp"   // mvie2d::inscribedEllipseTriangle(...)
#include "menn_quad_se.hpp"  // mvie2d::inscribedEllipseQuad(...)
#include "menn_pentagon.hpp" // mvie2d::inscribedEllipsePentagon(...)

namespace mvie2d
{

    // =============== 工具：面积/包含性/闭合性 ===============

    // 面积：π / sqrt(det(X))，等价于 E.area()
    inline double area_from_X(const Ellipse2 &E)
    {
        if (!E.valid)
            return 0.0;
        double detX = E.X.determinant();
        if (!(detX > 0.0))
            return 0.0;
        return M_PI / std::sqrt(detX);
    }

    // 半空间（ax+by+c<=0）类型，便于从顶点构造
    struct Halfspace
    {
        double a{0}, b{0}, c{0};
        Halfspace() = default;
        Halfspace(double A, double B, double C) : a(A), b(B), c(C) {}
        inline Eigen::Vector2d n() const { return {a, b}; }
        inline double eval(const Eigen::Vector2d &p) const { return a * p.x() + b * p.y() + c; }
        void normalize()
        {
            double l = std::hypot(a, b);
            if (l > 0)
            {
                a /= l;
                b /= l;
                c /= l;
            }
        }
    };

    // 采样检查椭圆是否在多边形内（以半空间集合表示）
    inline bool ellipse_inside_halfspaces_sampling(const Ellipse2 &E,
                                                   const std::vector<Halfspace> &H,
                                                   int samples = 128,
                                                   double tol = 1e-10)
    {
        if (!E.valid)
            return false;

        // X = SPD，解 X = L L^T
        Eigen::LLT<Eigen::Matrix2d> llt(E.X);
        if (llt.info() != Eigen::Success)
            return false;
        Eigen::Matrix2d L = llt.matrixL();

        for (int k = 0; k < samples; ++k)
        {
            double t = (2.0 * M_PI) * (double(k) / samples);
            Eigen::Vector2d u(std::cos(t), std::sin(t)); // 单位圆

            // 求 (x-c)^T X (x-c) = 1 上的点：
            // X = L L^T, 令 v = L^T (x-c) => ||v|| = 1
            // 取 v = u，则 (x-c) = (L^T)^{-1} u
            Eigen::Vector2d y = L.transpose().triangularView<Eigen::Upper>().solve(u);
            Eigen::Vector2d p = E.c_center + y;

            for (const auto &hs : H)
                if (hs.eval(p) > tol)
                    return false;
        }
        return true;
    }

    // 3/4/5 条切线是否构成“闭合”的简单多边形（按法向极角排）
    struct _Item
    {
        int idx;
        double ang;
        TangencyLine L;
    };

    inline bool subset_forms_closed_polygon(const std::vector<TangencyLine> &TL,
                                            const std::vector<int> &idxs,
                                            std::vector<TangencyLine> &ordered_out)
    {
        const int m = (int)idxs.size();
        if (m < 3 || m > 5)
            return false;

        std::vector<_Item> v;
        v.reserve(m);
        for (int id : idxs)
        {
            const auto &L = TL[id];
            double ang = std::atan2(L.b, L.a); // 法向极角
            v.push_back({id, ang, L});
        }
        std::sort(v.begin(), v.end(), [](const _Item &x, const _Item &y)
                  { return x.ang < y.ang; });

        auto intersect = [](const TangencyLine &L1, const TangencyLine &L2) -> std::optional<Eigen::Vector2d>
        {
            double det = L1.a * L2.b - L2.a * L1.b;
            if (std::fabs(det) < 1e-12)
                return std::nullopt;
            double x = (L2.b * (-L1.c) - L1.b * (-L2.c)) / det;
            double y = (L1.a * (-L2.c) - L2.a * (-L1.c)) / det;
            return Eigen::Vector2d(x, y);
        };

        std::vector<Eigen::Vector2d> verts;
        verts.reserve(m);
        for (int i = 0; i < m; ++i)
        {
            auto ip = intersect(v[i].L, v[(i + 1) % m].L);
            if (!ip)
                return false;
            verts.push_back(*ip);
        }

        // 顶点都要在所有半空间内（<=0）
        for (const auto &P : verts)
            for (const auto &it : v)
                if (it.L.a * P.x() + it.L.b * P.y() + it.L.c > 1e-8)
                    return false;

        ordered_out.clear();
        ordered_out.reserve(m);
        for (auto &it : v)
            ordered_out.push_back(it.L);

        return true;
    }

    // =============== 评价一个 3/4/5 子集：调用已有 MENN ===============

    struct SubsetEval
    {
        bool closed{false};
        bool ok{false};
        Ellipse2 E{};
        double area{0.0};
    };

    inline SubsetEval eval_subset_MENN(const std::vector<TangencyLine> &TL,
                                       const std::vector<int> &idxs,
                                       double eps = 1e-12)
    {
        SubsetEval out;
        std::vector<TangencyLine> ordered;
        if (!subset_forms_closed_polygon(TL, idxs, ordered))
        {
            out.closed = false;
            return out;
        }
        out.closed = true;

        // 3 条边：用三角形 Steiner 内接椭圆
        if (ordered.size() == 3)
        {
            Ellipse2 r = inscribedEllipseTriangle(ordered, eps); // 新接口：直接返回 Ellipse2
            if (r.valid)
            {
                out.ok = true;
                out.E = r;
                out.area = area_from_X(out.E);
            }
            return out;
        }

        // 4 条边：Horwitz 四边形最大面积椭圆
        if (ordered.size() == 4)
        {
            Ellipse2 r = inscribedEllipseQuad(ordered); // 新接口：直接返回 Ellipse2
            if (r.valid)
            {
                out.ok = true;
                out.E = r;
                out.area = area_from_X(out.E);
            }
            return out;
        }

        // 5 条边：五边形内接椭圆
        if (ordered.size() == 5)
        {
            Ellipse2 r = inscribedEllipsePentagon(ordered); // 新接口：直接返回 Ellipse2
            if (r.valid)
            {
                out.ok = true;
                out.E = r;
                out.area = area_from_X(out.E);
            }
            return out;
        }

        return out;
    }

    // =============== 顶点(CCW) → 切线(ax+by+c<=0 为内侧) ===============

    inline std::vector<TangencyLine> tangency_from_ccw_vertices(const std::vector<Eigen::Vector2d> &V,
                                                                double inward_eps = 0.0)
    {
        std::vector<TangencyLine> TL;
        const int n = (int)V.size();
        if (n < 3)
            return TL;

        TL.reserve(n);
        for (int i = 0; i < n; ++i)
        {
            Eigen::Vector2d p = V[i];
            Eigen::Vector2d q = V[(i + 1) % n];
            Eigen::Vector2d e = q - p;          // CCW 边方向
            Eigen::Vector2d nrm(e.y(), -e.x()); // 外法向
            Eigen::Vector2d inn = -nrm;         // 内法向
            double a = inn.x(), b = inn.y();
            double c = -(a * p.x() + b * p.y());
            double l = std::hypot(a, b);
            if (l > 0)
            {
                a /= l;
                b /= l;
                c /= l;
            }
            if (inward_eps != 0.0)
                c -= inward_eps;
            TL.emplace_back(a, b, c);
        }
        return TL;
    }

    // =============== 选优入口（半空间/切线集合） ===============

    struct AnyMVIEOptions
    {
        int random_trials = 2;    // 多次随机顺序
        int inside_samples = 128; // 包含性采样点
        double eps = 1e-12;
        unsigned seed = 114514;
    };

    // 主入口一：从切线集合（ax+by+c<=0）求任意凸多边形的“最大内接椭圆”
    inline Ellipse2 computeMaxArea_fromAnyConvexPolygon(const std::vector<TangencyLine> &TL_in,
                                                        const AnyMVIEOptions &opt = {})
    {
        Ellipse2 best;
        best.valid = false;
        if (TL_in.size() < 3)
            return best;

        // 小规范化一遍（不改变几何意义，仅改善数值）
        std::vector<TangencyLine> TL = TL_in;
        for (auto &L : TL)
        {
            double l = std::hypot(L.a, L.b);
            if (l > 0)
            {
                L.a /= l;
                L.b /= l;
                L.c /= l;
            }
        }

        std::mt19937 rng(opt.seed);
        std::vector<int> ids(TL.size());
        std::iota(ids.begin(), ids.end(), 0);

        auto build_halfspaces_all = [&]()
        {
            std::vector<Halfspace> H;
            H.reserve(TL.size());
            for (auto &t : TL)
            {
                Halfspace h{t.a, t.b, t.c};
                h.normalize();
                H.push_back(h);
            }
            return H;
        };

        auto try_subsets = [&](const std::vector<int> &order)
        {
            const int N = (int)order.size();
            std::vector<int> pick;

            // 全局 halfspace 集合（全集约束）
            auto H_all = build_halfspaces_all();

            // 3 组合
            for (int i = 0; i < N; ++i)
                for (int j = i + 1; j < N; ++j)
                    for (int k = j + 1; k < N; ++k)
                    {
                        pick = {order[i], order[j], order[k]};
                        auto ev = eval_subset_MENN(TL, pick, opt.eps);
                        if (!ev.closed || !ev.ok)
                            continue;
                        if (!ellipse_inside_halfspaces_sampling(ev.E, H_all, opt.inside_samples))
                            continue;
                        if (!best.valid || ev.area > area_from_X(best))
                            best = ev.E;
                    }

            // 4 组合
            for (int i = 0; i < N; ++i)
                for (int j = i + 1; j < N; ++j)
                    for (int k = j + 1; k < N; ++k)
                        for (int l = k + 1; l < N; ++l)
                        {
                            pick = {order[i], order[j], order[k], order[l]};
                            auto ev = eval_subset_MENN(TL, pick, opt.eps);
                            if (!ev.closed || !ev.ok)
                                continue;
                            if (!ellipse_inside_halfspaces_sampling(ev.E, H_all, opt.inside_samples))
                                continue;
                            if (!best.valid || ev.area > area_from_X(best))
                                best = ev.E;
                        }

            // 5 组合
            for (int i = 0; i < N; ++i)
                for (int j = i + 1; j < N; ++j)
                    for (int k = j + 1; k < N; ++k)
                        for (int l = k + 1; l < N; ++l)
                            for (int m = l + 1; m < N; ++m)
                            {
                                pick = {order[i], order[j], order[k], order[l], order[m]};
                                auto ev = eval_subset_MENN(TL, pick, opt.eps);
                                if (!ev.closed || !ev.ok)
                                    continue;
                                if (!ellipse_inside_halfspaces_sampling(ev.E, H_all, opt.inside_samples))
                                    continue;
                                if (!best.valid || ev.area > area_from_X(best))
                                    best = ev.E;
                            }
        };

        for (int t = 0; t < std::max(1, opt.random_trials); ++t)
        {
            std::shuffle(ids.begin(), ids.end(), rng);
            try_subsets(ids);
        }
        return best;
    }

    // mvie2d_any_polygon.hpp 末尾，namespace mvie2d 内

    // 模板入口：可以吃任意 TangencyLike，只要有 .a/.b/.c 成员
    template <class TangencyLike>
    inline Ellipse2 computeMaxArea_fromAnyConvexPolygon(
        const std::vector<TangencyLike> &TL_in,
        const AnyMVIEOptions &opt = {})
    {
        // 1) 先把用户的 TangencyLike（比如你的 Tangencyline）
        //    转成库内部用的精简 TangencyLine
        std::vector<TangencyLine> TL;
        TL.reserve(TL_in.size());
        for (const auto &t : TL_in)
        {
            TL.emplace_back(t.a, t.b, t.c);
        }

        // 2) 调用原来的实现（只操作 mvie2d::TangencyLine）
        return computeMaxArea_fromAnyConvexPolygon(TL, opt);
    }

    // 主入口二：从 CCW 顶点集合求 MVIE（内部自动转切线）
    inline Ellipse2 computeMaxArea_fromCCWVertices(const std::vector<Eigen::Vector2d> &V,
                                                   const AnyMVIEOptions &opt = {},
                                                   double inward_eps = 0.0)
    {
        auto TL = tangency_from_ccw_vertices(V, inward_eps);
        return computeMaxArea_fromAnyConvexPolygon(TL, opt);
    }

} // namespace mvie2d
