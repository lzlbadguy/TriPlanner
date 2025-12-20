#ifndef INTERSECT_COLLISION_H
#define INTERSECT_COLLISION_H

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <limits>
#include "utils.hpp"

typedef Eigen::Vector2d Vec2;

class node_with_expansion
{
public:
    double x;
    double y;
    double theta;
    double delta_s = 0.05;
    double L_limit = 6.0;
    std::vector<double> expanded_length = {3.88, 0.96, 1.12, 0.96};
    std::vector<int> expand_index = {0, 1, 2, 3};
    std::vector<int> erase_flag = {0, 0, 0, 0};
};

class Polygon_simple
{
public:
    std::vector<double> x;
    std::vector<double> y;
};

class EdgeIntersector
{
public:
    /**
     * Calculate the intersection of two groups of edges.
     *
     * @param e1 Coordinates of m edges in shape (m, 2, 2).
     * @param e2 Coordinates of n edges in shape (n, 2, 2).
     * @return Intersections: A matrix of shape (m, n, 2). If edges do not intersect, the corresponding value is set to infinity.
     */

    static void vertices_transfer(const node_with_expansion &node, Eigen::Matrix<double, 2, 4> &vehicle)
    {
        vehicle << node.expanded_length[0], node.expanded_length[0], -node.expanded_length[2], -node.expanded_length[2], -node.expanded_length[3], node.expanded_length[1], node.expanded_length[1], -node.expanded_length[3];
        Eigen::Matrix2d rot;
        rot << cos(node.theta), -sin(node.theta), sin(node.theta), cos(node.theta);
        vehicle = rot * vehicle;
        vehicle += Vec2(node.x, node.y).replicate(1, 4);
    }

    static void vehicle_coordinate_transfer_L_to_G(const node_with_expansion &node, Eigen::Matrix<double, 2, 4> &vehicle, const utils::VehicleConfig &vc)
    {   
        vehicle << vc.RF, vc.RF, -vc.RB, -vc.RB, -vc.W/2, vc.W/2, vc.W/2, -vc.W/2;
        Eigen::Matrix2d rot;
        rot << cos(node.theta), -sin(node.theta), sin(node.theta), cos(node.theta);
        vehicle = rot * vehicle;
        vehicle += Vec2(node.x, node.y).replicate(1, 4);
    }

    static bool collisioncheck(const Polygon_simple &poly1, const Polygon_simple &poly2)
    {
        std::vector<Eigen::Matrix2d> e1 = getedges(poly1.x, poly1.y);
        std::vector<Eigen::Matrix2d> e2 = getedges(poly2.x, poly2.y);
        return isintersect(e1, e2);
    }

    static std::vector<Eigen::Matrix2d> getedges(const std::vector<double> &x, const std::vector<double> &y)
    {
        std::vector<Eigen::Matrix2d> edges;
        for (size_t i = 0; i < x.size() - 1; i++)
        {
            Eigen::Matrix2d edge;
            edge << x[i], y[i], x[i + 1], y[i + 1];
            edges.push_back(edge);
        }
        return edges;
    }

    static bool isintersect(
        const std::vector<Eigen::Matrix2d> &e1,
        const std::vector<Eigen::Matrix2d> &e2)
    {
        size_t m = e1.size();
        size_t n = e2.size();

        // Initialize the result matrix with infinity
        // std::vector<std::vector<Eigen::Vector2d>> intersections(
        //     m, std::vector<Eigen::Vector2d>(n, Eigen::Vector2d::Constant(std::numeric_limits<double>::infinity()))
        // );
        bool iscollision = false;
        double tolerrance = 1e-8;

        for (size_t i = 0; i < m; ++i)
        {
            // Extract edge 1 parameters
            double x1s1 = e1[i](0, 0), y1s1 = e1[i](0, 1);
            double x2s1 = e1[i](1, 0), y2s1 = e1[i](1, 1);
            double a = y2s1 - y1s1;
            double b = x1s1 - x2s1;
            double c = y1s1 * x2s1 - x1s1 * y2s1;

            for (size_t j = 0; j < n; ++j)
            {
                // Extract edge 2 parameters
                double x1s2 = e2[j](0, 0), y1s2 = e2[j](0, 1);
                double x2s2 = e2[j](1, 0), y2s2 = e2[j](1, 1);
                double d = y2s2 - y1s2;
                double e = x1s2 - x2s2;
                double f = y1s2 * x2s2 - x1s2 * y2s2;

                // Compute determinant
                double det = a * e - b * d;
                if (std::abs(det) < tolerrance)
                {
                    // Lines are parallel
                    continue;
                }

                // Calculate intersection point
                double raw_x = (b * f - c * e) / det;
                double raw_y = (c * d - a * f) / det;

                // Check if the intersection point is within the bounds of both edges
                if (raw_x > std::max(x1s1, x2s1) + tolerrance || raw_x < std::min(x1s1, x2s1) - tolerrance ||
                    raw_y > std::max(y1s1, y2s1) + tolerrance || raw_y < std::min(y1s1, y2s1) - tolerrance ||
                    raw_x > std::max(x1s2, x2s2) + tolerrance || raw_x < std::min(x1s2, x2s2) - tolerrance ||
                    raw_y > std::max(y1s2, y2s2) + tolerrance || raw_y < std::min(y1s2, y2s2) - tolerrance)
                {
                    continue;
                }

                iscollision = true;
            }
        }

        return iscollision;
    }

    static void expandPolygon(node_with_expansion &node, const std::vector<Polygon_simple> &obstacles)
    {
        auto it = node.expand_index.begin();
        while (!node.expand_index.empty())
        {
            int direction = *it;
            auto last_expanded_length = node.expanded_length;
            node.expanded_length[direction] = node.expanded_length[direction] + node.delta_s;
            if (node.expanded_length[direction] >= node.L_limit)
            {
                it = node.expand_index.erase(it);
                node.expanded_length = last_expanded_length;
                if (it == node.expand_index.end())
                {
                    it = node.expand_index.begin();
                }
            }
            else
            {
                for (int j = 0; j < obstacles.size(); j++)
                {
                    Eigen::Matrix<double, 2, 4> vehicle;
                    vertices_transfer(node, vehicle);
                    Polygon_simple poly_box;
                    poly_box.x = {vehicle(0, 0), vehicle(0, 1), vehicle(0, 2), vehicle(0, 3), vehicle(0, 0)};
                    poly_box.y = {vehicle(1, 0), vehicle(1, 1), vehicle(1, 2), vehicle(1, 3), vehicle(1, 0)};
                    if (collisioncheck(poly_box, obstacles[j]))
                    {
                        it = node.expand_index.erase(it);
                        if (it == node.expand_index.end())
                        {
                            it = node.expand_index.begin();
                        }
                        node.expanded_length = last_expanded_length;
                        break;
                    }

                    if (j == obstacles.size() - 1)
                    {
                        ++it;
                        if (it == node.expand_index.end())
                        {
                            it = node.expand_index.begin();
                        }
                    }
                }
            }
        }
    }
};

#endif // INTERSECT_COLLISION_H