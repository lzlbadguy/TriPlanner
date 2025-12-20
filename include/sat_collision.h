#ifndef SAT_COLLISION_H
#define SAT_COLLISION_H

#include <Eigen/Dense>
#include <vector>
#include <limits>
#include <cmath>

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
    // std::vector<double> expanded_length = {0.8, 0.8, 0.8, 0.8};
    std::vector<int> expand_index = {0, 1, 2, 3};
    std::vector<int> erase_flag = {0, 0, 0, 0};
};

class Polygon
{
public:
    std::vector<Vec2> vertices;

    Polygon(const std::vector<Vec2> &verts) : vertices(verts) {}

    // 获取边列表
    Eigen::MatrixXd getEdges() const
    {
        Eigen::MatrixXd edges(2, vertices.size());
        for (size_t i = 0; i < vertices.size(); ++i)
        {
            Vec2 edge = vertices[(i + 1) % vertices.size()] - vertices[i];
            edges.col(i) = edge;
        }
        return edges;
    }
};

class SATCollision
{
public:
    // 投影多边形到轴
    static void projectPolygon(const Polygon &polygon, const Vec2 &axis, double &min, double &max)
    {
        min = std::numeric_limits<double>::infinity();
        max = -std::numeric_limits<double>::infinity();

        for (const Vec2 &vertex : polygon.vertices)
        {
            double projection = vertex.dot(axis);
            if (projection < min)
                min = projection;
            if (projection > max)
                max = projection;
        }
    }

    // 检查投影是否重叠
    static bool overlapOnAxis(const Polygon &poly1, const Polygon &poly2, const Vec2 &axis)
    {
        double min1, max1, min2, max2;
        projectPolygon(poly1, axis, min1, max1);
        projectPolygon(poly2, axis, min2, max2);

        return !(max1 < min2 || max2 < min1); // 如果没有重叠，则返回 false
    }

    // SAT 碰撞检测主函数
    static bool isColliding(const Polygon &poly1, const Polygon &poly2)
    {
        Eigen::MatrixXd edges1 = poly1.getEdges();
        Eigen::MatrixXd edges2 = poly2.getEdges();

        // 检查所有边的法线（即分离轴）
        for (size_t i = 0; i < edges1.cols(); ++i)
        {
            Vec2 edge = edges1.col(i);
            Vec2 axis(-edge.y(), edge.x()); // 获取边的法线
            axis.normalize();
            if (!overlapOnAxis(poly1, poly2, axis))
                return false; // 如果某个轴上没有重叠
        }

        for (size_t i = 0; i < edges2.cols(); ++i)
        {
            Vec2 edge = edges2.col(i);
            Vec2 axis(-edge.y(), edge.x()); // 获取边的法线
            axis.normalize();
            if (!overlapOnAxis(poly1, poly2, axis))
                return false; // 如果某个轴上没有重叠
        }

        return true; // 所有轴上都有重叠，说明碰撞
    }

    static void vertices_transfer(const node_with_expansion &node, Eigen::Matrix<double, 2, 4> &vehicle)
    {
        vehicle << node.expanded_length[0], node.expanded_length[0], -node.expanded_length[2], -node.expanded_length[2], -node.expanded_length[3], node.expanded_length[1], node.expanded_length[1], -node.expanded_length[3];
        Eigen::Matrix2d rot;
        rot << cos(node.theta), -sin(node.theta), sin(node.theta), cos(node.theta);
        vehicle = rot * vehicle;
        vehicle += Vec2(node.x, node.y).replicate(1, 4);
    }

    static void expandPolygon(node_with_expansion &node, const std::vector<Polygon> &obstacles)
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
                    Polygon poly_box({Vec2(vehicle(0, 0), vehicle(1, 0)),
                                      Vec2(vehicle(0, 1), vehicle(1, 1)),
                                      Vec2(vehicle(0, 2), vehicle(1, 2)),
                                      Vec2(vehicle(0, 3), vehicle(1, 3))});
                    // auto x_coord = node.x;
                    // auto y_coord = node.y;
                    // auto el = node.expanded_length;
                    // Polygon poly_box({Vec2(x_coord + el[0], y_coord - el[3]),
                    //                   Vec2(x_coord + el[0], y_coord + el[1]),
                    //                   Vec2(x_coord - el[2], y_coord + el[1]),
                    //                   Vec2(x_coord - el[2], y_coord - el[3])});
                    if (isColliding(poly_box, obstacles[j]))
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

#endif // SAT_COLLISION_H
