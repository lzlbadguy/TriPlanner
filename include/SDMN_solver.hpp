#ifndef SDMN_SOLVER_HPP
#define SDMN_SOLVER_HPP

#include <Eigen/Dense>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <cassert>

class SDMNSolver
{
public:
    using Vec2 = Eigen::Vector2d;
    using Mat2 = Eigen::Matrix2d;

    struct Halfspace
    {
        Vec2 a;
        double b;
    };

    struct Halfspace_1D
    {
        double a;
        double b;
    };

    static Vec2 findTangencyPoint(const std::vector<Vec2> &obs_pts, const std::vector<Vec2> &veh_pts)
    {
        std::vector<Halfspace> constraints;
        for (const auto &pt : veh_pts)
        {
            constraints.push_back({pt, 1.0});
        }

        for (const auto &pt : obs_pts)
        {
            constraints.push_back({-pt, -1.0});
        }

        Vec2 beta = SDMN(constraints);

        beta /= beta.squaredNorm();

        return beta;
    }

private:
    static bool isViolated(const Vec2 &y, const Halfspace &h)
    {
        constexpr double kTolerance = 1e-9;
        return h.a.dot(y) > h.b + kTolerance;
    }

    static double SolveOneDimMinNorm(const std::vector<Halfspace_1D> &H_reduced)
    {
        double y_min = -std::numeric_limits<double>::infinity();
        double y_max = std::numeric_limits<double>::infinity();

        for (const auto &h : H_reduced)
        {
            double a = h.a;
            double b = h.b;

            if (std::abs(a) < 1e-9)
            {
                if (b < -1e-9)
                    return std::numeric_limits<double>::infinity(); // infeasible
                continue;
            }

            double bound = b / a;
            if (a > 0)
            {
                y_max = std::min(y_max, bound);
            }
            else
            {
                y_min = std::max(y_min, bound);
            }
        }

        if (y_min > y_max)
        {
            return std::numeric_limits<double>::infinity(); // infeasible
        }

        // min y'^2
        if (y_min <= 0 && 0 <= y_max)
        {
            return 0.0;
        }
        else if (std::abs(y_min) < std::abs(y_max))
        {
            return y_min;
        }
        else
        {
            return y_max;
        }
    }

    static void HouseholderProj(const std::vector<Halfspace> &I, const Halfspace &h,
                                Vec2 &M, Vec2 &v, std::vector<Halfspace_1D> &H_reduced)
    {
        const Vec2 &a = h.a;
        double b = h.b;

        v = (b / a.squaredNorm()) * a;

        int j;
        v.cwiseAbs().maxCoeff(&j);
        Vec2 u = v;
        u[j] += (v[j] >= 0 ? 1 : -1) * v.norm();

        u.normalize();

        Mat2 H = Mat2::Identity() - 2.0 * u * u.transpose();

        M = H.transpose().col(1 - j);

        H_reduced.clear();
        for (const auto &hi : I)
        {
            Halfspace_1D h_new;
            h_new.a = M.dot(hi.a);
            h_new.b = hi.b - hi.a.dot(v);
            H_reduced.push_back(h_new); // 修复变量名：h_reduced → H_reduced
        }
    }

    static Vec2 SDMN(const std::vector<Halfspace> &H_E)
    {
        // Vec2 y = Vec2::Ones(); // instead of Zero
        Vec2 y = Vec2::Zero();
        std::vector<Halfspace> I;
        std::vector<Halfspace> H = H_E;
        std::shuffle(H.begin(), H.end(), std::mt19937{std::random_device{}()});

        for (const auto &h : H)
        {
            if (!isViolated(y, h))
            {
                I.push_back(h);
                continue;
            }

            Vec2 M, v;
            std::vector<Halfspace_1D> H_reduced;
            HouseholderProj(I, h, M, v, H_reduced);

            double y_prime = 0.0;
            if (!H_reduced.empty())
            {
                y_prime = SolveOneDimMinNorm(H_reduced);
            }

            Vec2 y_full = M * y_prime + v;
            y = y_full;
            I.push_back(h);
        }

        return y;
    }
};

#endif // SDMN_SOLVER_HPP
