#include <iostream>
#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <cmath>
#include <numeric>
#include "matplotlibcpp.h"
#include "utils.hpp"
#include "SDMN_solver.hpp"
#include "mvie2d_types.hpp"
#include "mvie2d_any_polygon.hpp"
using namespace Eigen;
typedef Eigen::Vector2d Vec2;
namespace plt = matplotlibcpp;

class simple_vehicle_with_position
{
public:
    double LF, LR, W, x, y, theta;
    std::vector<Vec2> corner_points;
    std::vector<double> corner_points_x, corner_points_y;
    Matrix2d R;
    simple_vehicle_with_position(double LF_, double LR_, double W_, double x_, double y_, double theta_)
        : LF(LF_), LR(LR_), W(W_), x(x_), y(y_), theta(theta_)
    {
        R << cos(theta), -sin(theta),
            sin(theta), cos(theta);
        corner_points.push_back(Vec2(LF, -W / 2));
        corner_points.push_back(Vec2(LF, W / 2));
        corner_points.push_back(Vec2(-LR, W / 2));
        corner_points.push_back(Vec2(-LR, -W / 2));
        for (int i = 0; i < 4; i++)
        {
            corner_points[i] = R * corner_points[i] + Vec2(x, y);
            corner_points_x.push_back(corner_points[i][0]);
            corner_points_y.push_back(corner_points[i][1]);
        }
    }
};

class Ellipse_with_parameters
{
public:
    double a, b, theta_r;
    Vec2 center;
    Matrix2d R, R_inv;

    // **默认构造函数**
    Ellipse_with_parameters()
        : a(1.0), b(1.0), theta_r(0.0), center(Vec2::Zero()), R(Matrix2d::Identity()), R_inv(Matrix2d::Identity())
    {
        // 初始化 R 和 R_inv
        R << cos(theta_r), -sin(theta_r),
            sin(theta_r), cos(theta_r);

        R_inv << cos(theta_r), sin(theta_r),
            -sin(theta_r), cos(theta_r);
    }
    // 构造函数
    Ellipse_with_parameters(double a_, double b_, double theta_r_, Vec2 center_)
        : a(a_), b(b_), theta_r(theta_r_), center(center_)
    {
        R << cos(theta_r), -sin(theta_r),
            sin(theta_r), cos(theta_r);

        R_inv << cos(theta_r), sin(theta_r),
            -sin(theta_r), cos(theta_r);
    }

    Ellipse_with_parameters &operator=(const Ellipse_with_parameters &other)
    {
        if (this != &other)
        {
            a = other.a;
            b = other.b;
            theta_r = other.theta_r;
            center = other.center;
            R = other.R;
            R_inv = other.R_inv;
        }
        return *this;
    }

    Vec2 transfer2unit_circle(const Vec2 &p) const
    {
        Vec2 p_rotated = R_inv * (p - center);
        Vec2 p_unit_circle(p_rotated(0) / a, p_rotated(1) / b);
        return p_unit_circle;
    }

    Vec2 unit_circle2Cartesian(const Vec2 &p) const
    {
        Vec2 p_reverse1(a * p(0), b * p(1));
        Vec2 p_cartesian = R * p_reverse1 + center;
        return p_cartesian;
    }
};

class Obstacle
{
public:
    std::vector<double> x, y;
};

bool check_opt_necessity(const Vec2 &tangency_point, const std::vector<double> &trasfer_obs_x1, const std::vector<double> &trasfer_obs_y1, const std::vector<double> &tf_vehicle_x, const std::vector<double> &tf_vehicle_y)
{
    auto a = tangency_point(0);
    auto b = tangency_point(1);
    auto c = -(a * a + b * b);
    bool Need_optimization = false;
    for (int i = 0; i < tf_vehicle_x.size(); i++)
    {
        auto f_vehicle = a * tf_vehicle_x[i] + b * tf_vehicle_y[i] + c;
        if (f_vehicle > 0.0)
        {
            Need_optimization = true;
            break;
        }
    }

    if (!Need_optimization)
    {
        for (int i = 0; i < trasfer_obs_x1.size(); i++)
        {
            auto f_obs = a * trasfer_obs_x1[i] + b * trasfer_obs_y1[i] + c;
            if (f_obs < 0.0)
            {
                Need_optimization = true;
                break;
            }
        }
    }

    return Need_optimization;
}

class Tangencyline
{
public:
    double a, b, c;
    // bool is_valid;
    Vec2 tangency_point_cartesian;
    Vec2 tangent_dir_ellipse;
    Vec2 ellipse_center_to_tangency_point;
    double tangency_point_to_ellipse_center_distance;
    std::vector<Vec2> tf_veh_pts, tf_obs_pts;
    Ellipse_with_parameters ellipse;
    Obstacle Obs;
    simple_vehicle_with_position vehicle;
    // 比较两个路径对象的代价
    bool operator<(const Tangencyline &other) const { return tangency_point_to_ellipse_center_distance < other.tangency_point_to_ellipse_center_distance; }
    // 构造函数（直接计算切线）
    Tangencyline(const Ellipse_with_parameters &ellipse_, const Obstacle &obs_, const simple_vehicle_with_position &vehicle_)
        : ellipse(ellipse_), Obs(obs_), vehicle(vehicle_), a(0), b(0), c(0)
    {
        // 1. 变换障碍物到单位圆坐标系
        std::vector<double> trasfer_obs_x1, trasfer_obs_y1, tf_vehicle_x, tf_vehicle_y;
        for (size_t k = 0; k < Obs.x.size(); k++)
        {
            Vec2 p_obs(Obs.x[k], Obs.y[k]);
            Vec2 obs1_unit_circle = ellipse.transfer2unit_circle(p_obs);
            trasfer_obs_x1.push_back(obs1_unit_circle(0));
            trasfer_obs_y1.push_back(obs1_unit_circle(1));
            tf_obs_pts.push_back(obs1_unit_circle);
        }

        for (size_t k = 0; k < vehicle.corner_points_x.size(); k++)
        {
            Vec2 p_vehicle(vehicle.corner_points_x[k], vehicle.corner_points_y[k]);
            Vec2 vehicle1_unit_circle = ellipse.transfer2unit_circle(p_vehicle);
            tf_vehicle_x.push_back(vehicle1_unit_circle(0));
            tf_vehicle_y.push_back(vehicle1_unit_circle(1));
            tf_veh_pts.push_back(vehicle1_unit_circle);
        }

        // 2. 计算单位圆坐标系下的切点
        Vec2 tangency_point = SDMNSolver::findTangencyPoint(tf_obs_pts, tf_veh_pts);
        // 3. 转换切点回原椭圆坐标系
        tangency_point_cartesian = ellipse.unit_circle2Cartesian(tangency_point);

        // 4. 计算切线方向向量（单位圆）
        Vec2 tangent_direction(-tangency_point(1), tangency_point(0)); // 切线方向 1
        tangent_dir_ellipse = ellipse.R * Vec2(ellipse.a * tangent_direction(0), ellipse.b * tangent_direction(1));
        ellipse_center_to_tangency_point = ellipse.center - tangency_point_cartesian;
        // tangency_point_to_ellipse_center_distance = ellipse_center_to_tangency_point.norm();
        tangency_point_to_ellipse_center_distance = tangency_point.norm();
        // 5. 计算切线方程 ax + by + c = 0
        Vec2 dir = tangent_dir_ellipse;
        double dx = dir(0), dy = dir(1);
        Vec2 P = tangency_point_cartesian;

        // 法向量为 (dy, -dx)，直线方程 dy*(x - P.x()) - dx*(y - P.y()) = 0
        a = dy;
        b = -dx;
        c = dx * P.y() - dy * P.x();
        // 标准化法向量朝向，使椭圆中心在“法向量负方向”（内侧）
        double f_center = a * ellipse.center.x() + b * ellipse.center.y() + c;
        if (f_center > 0) // 椭圆中心在法向量“正”方向 → 说明方向反了，需要flip
        {
            a *= -1;
            b *= -1;
            c *= -1;
        }
    }

    bool inside_half_plain(const double &x, const double &y, const simple_vehicle_with_position &vehicle)
    {
        double f = a * x + b * y + c;
        double f_vehicle = a * vehicle.x + b * vehicle.y + c;

        return (f * f_vehicle < 0);
    }
};

bool within_all_Tangencylines(const std::vector<Tangencyline> &Tangency_lines, const Vec2 &intersection)
{
    for (const auto &line : Tangency_lines)
    {
        double f = line.a * intersection.x() + line.b * intersection.y() + line.c;
        if (f > 1e-6) // >0 表示在切线外侧，剔除；允许一定数值容差
            return false;
    }
    return true;
}

// 辅助函数：计算两条切线的交点
bool computeIntersection(const Tangencyline &line1, const Tangencyline &line2, Vec2 &intersection)
{
    double a1 = line1.a, b1 = line1.b, c1 = line1.c;
    double a2 = line2.a, b2 = line2.b, c2 = line2.c;

    double det = a1 * b2 - a2 * b1;
    if (std::abs(det) < 1e-6)
    {
        return false; // 平行或重合
    }

    double x = (b1 * c2 - b2 * c1) / det;
    double y = (a2 * c1 - a1 * c2) / det;
    intersection = Vec2(x, y);
    return true;
}

std::vector<Tangencyline> FilterOut_Invalid_Tangencyline(std::vector<Tangencyline> &Tangency_lines, const simple_vehicle_with_position &vehicle)
{
    std::vector<Tangencyline> valid_Tangencylines;
    std::sort(Tangency_lines.begin(), Tangency_lines.end());
    std::vector<int> invalid_indexs;
    for (int i = 0; i < Tangency_lines.size(); i++)
    {
        auto it = std::find(invalid_indexs.begin(), invalid_indexs.end(), i);
        if (it != invalid_indexs.end())
        {
            continue;
        }
        else
        {
            for (int j = i + 1; j < Tangency_lines.size(); j++)
            {
                Obstacle check_obs = Tangency_lines[j].Obs;
                for (int k = 0; k < check_obs.x.size(); k++)
                {
                    if (!Tangency_lines[i].inside_half_plain(check_obs.x[k], check_obs.y[k], vehicle))
                    {
                        break;
                    }
                    else
                    {
                        if (k == check_obs.x.size() - 1)
                        {
                            auto it = std::find(invalid_indexs.begin(), invalid_indexs.end(), j);
                            if (it == invalid_indexs.end())
                            {
                                invalid_indexs.push_back(j);
                                // std::cout << "Tangencyline " << i << " removed the invalid Tangencyline " << j << std::endl;
                            }

                            break;
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < Tangency_lines.size(); i++)
    {
        auto it = std::find(invalid_indexs.begin(), invalid_indexs.end(), i);
        if (it != invalid_indexs.end())
        {
            continue;
        }
        else
        {
            valid_Tangencylines.push_back(Tangency_lines[i]);
        }
    }

    return valid_Tangencylines;
}

void Generate_polygonVertices_basedon_Tangencyline(const Ellipse_with_parameters &ellipse, std::vector<Tangencyline> &Tangency_lines, std::vector<Vec2> &polygon_vertices)
{
    std::vector<int> valid_indexs;
    std::vector<Tangencyline> valid_Tangencylines;
    // 1. 按切点的极角排序
    std::sort(Tangency_lines.begin(), Tangency_lines.end(), [&ellipse](const Tangencyline &a, const Tangencyline &b)
              {
      Vec2 rel_a = a.tangency_point_cartesian - ellipse.center;
      double angle_a = atan2(rel_a.y(), rel_a.x());
      Vec2 rel_b = b.tangency_point_cartesian - ellipse.center;
      double angle_b = atan2(rel_b.y(), rel_b.x());
      return angle_a < angle_b; });

    for (size_t i = 0; i < Tangency_lines.size(); ++i)
    {
        Tangencyline current = Tangency_lines[i];
        for (int j = i + 1; j < Tangency_lines.size(); ++j)
        {

            Tangencyline next = Tangency_lines[j];
            Vec2 intersection;
            if (computeIntersection(current, next, intersection))
            {
                if (within_all_Tangencylines(Tangency_lines, intersection))
                {
                    polygon_vertices.push_back(intersection);
                    auto it = std::find(valid_indexs.begin(), valid_indexs.end(), i);
                    if (it == valid_indexs.end())
                    {
                        valid_indexs.push_back(i);
                        valid_Tangencylines.push_back(current);
                    }
                    auto it2 = std::find(valid_indexs.begin(), valid_indexs.end(), j);
                    if (it2 == valid_indexs.end())
                    {
                        valid_indexs.push_back(j);
                        valid_Tangencylines.push_back(next);
                    }
                }
            }
        }
    }

    Tangency_lines = valid_Tangencylines;
}

std::vector<Tangencyline> Generate_Edges_for_single_safety_block(const Ellipse_with_parameters &ellipse, const std::vector<Obstacle> &obstacles, const simple_vehicle_with_position &vehicle)
{
    std::vector<Tangencyline> Tangency_lines;
    std::vector<Obstacle> obstacles_within_range;
    // 创建 Tangencyline 对象（自动计算切线）
    for (int i = 0; i < obstacles.size(); i++)
    {
        if (i == obstacles.size() - 1 || i == obstacles.size() - 2)
        {
            obstacles_within_range.push_back(obstacles[i]);
            continue;
        }

        for (int j = 0; j < obstacles[i].x.size(); j++)
        {
            double dist = (Vec2(obstacles[i].x[j], obstacles[i].y[j]) - Vec2(vehicle.x, vehicle.y)).norm();
            if (dist < 30.0)
            {
                obstacles_within_range.push_back(obstacles[i]);
                break;
            }
        }
    }

    for (int i = 0; i < obstacles_within_range.size(); i++)
    {
        Tangencyline tangent_line(ellipse, obstacles_within_range[i], vehicle);
        Tangency_lines.push_back(tangent_line);
    }

    std::vector<Tangencyline> valid_Tangencylines = FilterOut_Invalid_Tangencyline(Tangency_lines, vehicle);
    return valid_Tangencylines;
}

Ellipse_with_parameters GenerateMaxInscribedEllipse(const std::vector<Tangencyline> &Tangency_lines)
{
    // 求解最大内接椭圆
    mvie2d::Ellipse2 E = mvie2d::computeMaxArea_fromAnyConvexPolygon(Tangency_lines);
    Ellipse_with_parameters Max_Inscribed_Ellipse(E.major_axis, E.minor_axis, E.c_theta, E.c_center);
    return Max_Inscribed_Ellipse;
}

void Generate_points_on_Ellipse(const Ellipse_with_parameters &ellipse, const int &num_points, std::vector<double> &x_vals, std::vector<double> &y_vals)
{
    for (int i = 0; i <= num_points; ++i)
    {
        double theta = 2.0 * M_PI * i / num_points; // 角度参数
        Vec2 p(cos(theta), sin(theta));             // 标准椭圆上的点
        Vec2 p_ellipse = ellipse.unit_circle2Cartesian(p);
        x_vals.push_back(p_ellipse(0));
        y_vals.push_back(p_ellipse(1));
    }
}

void Iterate_safety_block(const Ellipse_with_parameters &ellipse, const std::vector<Obstacle> &obstacles, const int &Max_iteration_num, const double &stop_iteration_criterion, std::vector<Tangencyline> &Tangency_lines, Ellipse_with_parameters &Max_Inscribed_Ellipse, const simple_vehicle_with_position &vehicle)
{
    // 创建 Tangencyline 对象（自动计算切线）
    for (int i = 0; i < Max_iteration_num; i++)
    {
        if (i == 0)
        {
            Tangency_lines = Generate_Edges_for_single_safety_block(ellipse, obstacles, vehicle);
            Max_Inscribed_Ellipse = GenerateMaxInscribedEllipse(Tangency_lines);
        }
        else
        {
            std::vector<Tangencyline> Tangency_lines_new;
            Tangency_lines_new = Generate_Edges_for_single_safety_block(Max_Inscribed_Ellipse, obstacles, vehicle);
            Ellipse_with_parameters Max_Inscribed_Ellipse_new = GenerateMaxInscribedEllipse(Tangency_lines_new);

            double proportion_of_expansion = (Max_Inscribed_Ellipse_new.a * Max_Inscribed_Ellipse_new.b - Max_Inscribed_Ellipse.a * Max_Inscribed_Ellipse.b) / (Max_Inscribed_Ellipse.a * Max_Inscribed_Ellipse.b);

            if (proportion_of_expansion < stop_iteration_criterion)
            {
                Max_Inscribed_Ellipse = Max_Inscribed_Ellipse_new;
                break;
            }
            else
            {
                Tangency_lines = Tangency_lines_new;
                Max_Inscribed_Ellipse = Max_Inscribed_Ellipse_new;
            }
        }
    }
}

Vec2 Convert_rear_axis_center_2_ellipse_center(const simple_vehicle_with_position &vehicle)
{
    Vec2 center;
    Vec2 rear_axis_to_center((vehicle.LF + vehicle.LR) / 2.0 - vehicle.LR, 0.0);
    Vec2 rear_axis(vehicle.x, vehicle.y);
    Matrix2d R;
    R << cos(vehicle.theta), -sin(vehicle.theta),
        sin(vehicle.theta), cos(vehicle.theta);
    center = R * rear_axis_to_center + rear_axis;
    return center;
}

class Parameters_for_optimization
{
public:
    int Max_iteration_num = 10;
    int num_points = 100;
    double stop_iteration_criterion = 1.0 / 100.0;
};

void Optimization_flow(const Ellipse_with_parameters &ellipse, const std::vector<Obstacle> &obstacles, const Parameters_for_optimization &params, std::vector<Tangencyline> &Tangency_lines, Ellipse_with_parameters &Max_Inscribed_Ellipse, const simple_vehicle_with_position &han, std::vector<Vec2> &polygon_vertices)
{
    Iterate_safety_block(ellipse, obstacles, params.Max_iteration_num, params.stop_iteration_criterion, Tangency_lines, Max_Inscribed_Ellipse, han);
    Generate_polygonVertices_basedon_Tangencyline(ellipse, Tangency_lines, polygon_vertices);
}

double polygonArea(const std::vector<double> &x, const std::vector<double> &y)
{
    // if (x.size() != y.size() || x.size() < 3) {
    //     throw std::invalid_argument("多边形至少需要3个顶点，且x和y大小要一致");
    // }
    const size_t n = x.size();
    double area = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        size_t j = (i + 1) % n;
        area += x[i] * y[j];
        area -= x[j] * y[i];
    }
    area = std::abs(area) * 0.5;
    return area;
}

void save_to_csv(const std::string &filename, const std::vector<Tangencyline> &Tangency_lines)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // 写入表头
    file << "a,b,c\n";

    // 写入每一行数据
    for (size_t i = 0; i < Tangency_lines.size(); i++)
    {
        file << Tangency_lines[i].a << "," << Tangency_lines[i].b << "," << Tangency_lines[i].c << "\n";
    }

    file.close();
    std::cout << "数据已写入到: " << filename << std::endl;
}

void save_to_csv_xy(const std::string &filename, const std::vector<double> &x, const std::vector<double> &y)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // 写入表头
    file << "x,y\n";

    // 写入每一行数据
    for (size_t i = 0; i < x.size(); i++)
    {
        file << x[i] << "," << y[i] << "\n";
    }

    file.close();
    std::cout << "数据已写入到: " << filename << std::endl;
}

void save_to_csv_vertices(const std::string &filename, const std::vector<Vec2> &polygon_vertices)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // 写入表头
    file << "x,y\n";

    // 写入每一行数据
    for (size_t i = 0; i < polygon_vertices.size() + 1; i++)
    {

        if (i == polygon_vertices.size())
        {
            file << polygon_vertices[0].x() << "," << polygon_vertices[0].y() << "\n";
        }
        else
        {
            file << polygon_vertices[i].x() << "," << polygon_vertices[i].y() << "\n";
        }
    }

    file.close();
}

int main()
{
    utils::file_manager fm;
    std::vector<double> ax, ay, ath, x_lattice, y_lattice, theta_lattice;
    Obstacle obs_vehicle, obs_vehicle1, obs_lower_bound, obs_upper_bound;
    std::vector<Obstacle> obstacles;
    // 基础路径
    std::string base_path = "../Data/Path_with_triangle_obstacles/";

    // 读取 Predicted_trajectory
    std::string predicted_trajectory_file = base_path + "reference_track.csv";
    std::string lattice_path_file = base_path + "lattice_path.csv";
    fm.read_csv_with_theta(predicted_trajectory_file, ax, ay, ath);
    fm.read_csv_with_theta(lattice_path_file, x_lattice, y_lattice, theta_lattice);
    // 读取 obs 数据
    for (int i = 0; i <= 19; ++i)
    {
        std::string obs_file = base_path + "obs_" + std::to_string(i) + ".csv";
        Obstacle obstacle;
        fm.read_csv(obs_file, obstacle.x, obstacle.y);
        obstacles.push_back(obstacle);
    }

    const double LF = 3.88;
    const double LR = 1.12;
    const double W = 1.92;
    double a = (LF + LR) / 2.0; // 长轴
    double b = W / 2.0;         // 短轴
    utils::TicToc Timer;

    simple_vehicle_with_position han_obs(LF, LR, W, ax[10], ay[10] - 1.5, ath[10]);
    simple_vehicle_with_position han_obs1(LF, LR, W, ax[20], ay[20] + 1.5, ath[10]);
    obs_vehicle.x = han_obs.corner_points_x;
    obs_vehicle.x.push_back(obs_vehicle.x[0]);
    obs_vehicle.y = han_obs.corner_points_y;
    obs_vehicle.y.push_back(obs_vehicle.y[0]);
    obs_vehicle1.x = han_obs1.corner_points_x;
    obs_vehicle1.x.push_back(obs_vehicle1.x[0]);
    obs_vehicle1.y = han_obs1.corner_points_y;
    obs_vehicle1.y.push_back(obs_vehicle1.y[0]);
    obs_lower_bound.x = {-5.0, -5.0};
    obs_lower_bound.y = {-5.0, 10.0};
    obs_upper_bound.x = {55.0, 55.0};
    obs_upper_bound.y = {5.0, 15.0};
    obstacles.push_back(obs_vehicle);
    obstacles.push_back(obs_vehicle1);
    obstacles.push_back(obs_lower_bound);
    obstacles.push_back(obs_upper_bound);

    Parameters_for_optimization params;
    std::vector<double> rec_cor_area, cor_area, cor_area_no_iteration;
    for (int i = 0; i < x_lattice.size(); i++)
    {
        plt::cla();
        // 1. 绘制障碍物
        simple_vehicle_with_position han(LF, LR, W, x_lattice[i], y_lattice[i], theta_lattice[i]);
        Vec2 center = Convert_rear_axis_center_2_ellipse_center(han);

        Ellipse_with_parameters ellipse(a, b, han.theta, center);

        std::vector<double> x_vals, y_vals;
        Generate_points_on_Ellipse(ellipse, params.num_points, x_vals, y_vals);

        std::vector<Tangencyline> Tangency_lines;
        Ellipse_with_parameters Max_Inscribed_Ellipse;
        double Calculate_time = 0.0;
        Timer.tic();
        std::vector<Vec2> polygon_vertices;
        Optimization_flow(ellipse, obstacles, params, Tangency_lines, Max_Inscribed_Ellipse, han, polygon_vertices);
        std::string Safety_Corridor_Path = "../Data/Path_with_triangle_obstacles/Safety_corridor/";
        std::string Safety_Corridor_Path_no_iteration = "../Data/Path_with_triangle_obstacles/Safety_corridor_no_iteration/";
        std::string filename = Safety_Corridor_Path + "Tangency_lines_" + std::to_string(i) + ".csv";
        std::string filename_no_iteration = Safety_Corridor_Path_no_iteration + "Tangency_lines_" + std::to_string(i) + ".csv";
        std::string filename_vertices = Safety_Corridor_Path + "Polygon_vertices_" + std::to_string(i) + ".csv";
        std::string filename_vertices_no_iteration = Safety_Corridor_Path_no_iteration + "Polygon_vertices_" + std::to_string(i) + ".csv";

        if (params.Max_iteration_num == 1)
        {
            save_to_csv(filename_no_iteration, Tangency_lines);
        }
        else
        {
            save_to_csv(filename, Tangency_lines);
        }

        Calculate_time = Timer.toc();
        std::cout << "Calculate_time: " << Calculate_time << "ms" << std::endl;
        if (Tangency_lines.size() == polygon_vertices.size())
        {
            std::cout << "Tangency_lines.size() == polygon_vertices.size()" << std::endl;
        }
        else
        {
            std::cout << "Tangency_lines.size() != polygon_vertices.size()" << std::endl;
        }

        std::sort(polygon_vertices.begin(), polygon_vertices.end(), [&ellipse](const Vec2 &a, const Vec2 &b)
                  {
        Vec2 rel_a = a - ellipse.center;
        double angle_a = atan2(rel_a.y(), rel_a.x());
        Vec2 rel_b = b - ellipse.center;
        double angle_b = atan2(rel_b.y(), rel_b.x());
        return angle_a < angle_b; });

        if (params.Max_iteration_num == 1)
        {
            save_to_csv_vertices(filename_vertices_no_iteration, polygon_vertices);
        }
        else
        {
            save_to_csv_vertices(filename_vertices, polygon_vertices);
        }

        // 3. 生成凸多边形的顶点数据
        std::vector<double> poly_x, poly_y;
        for (const auto &pt : polygon_vertices)
        {
            poly_x.push_back(pt(0));
            poly_y.push_back(pt(1));
        }
        if (!poly_x.empty())
        { // 闭合多边形
            poly_x.push_back(polygon_vertices[0](0));
            poly_y.push_back(polygon_vertices[0](1));
        }

        std::vector<double> c_x_vals, c_y_vals;
        Generate_points_on_Ellipse(Max_Inscribed_Ellipse, params.num_points, c_x_vals, c_y_vals);
        han.corner_points_x.push_back(han.corner_points_x[0]);
        han.corner_points_y.push_back(han.corner_points_y[0]);

        // for (int i = 0; i < Tangency_lines.size(); i++)
        // {
        //     std::vector<Vec2> points_on_tangencyline;
        //     std::vector<double> x_vals_tangencyline, y_vals_tangencyline;
        //     points_on_tangencyline.push_back(Tangency_lines[i].tangency_point_cartesian - 3 * Tangency_lines[i].tangent_dir_ellipse);
        //     points_on_tangencyline.push_back(Tangency_lines[i].tangency_point_cartesian + 3 * Tangency_lines[i].tangent_dir_ellipse);
        //     for (int j = 0; j < points_on_tangencyline.size(); j++)
        //     {
        //         x_vals_tangencyline.push_back(points_on_tangencyline[j].x());
        //         y_vals_tangencyline.push_back(points_on_tangencyline[j].y());
        //     }
        //     plt::plot(x_vals_tangencyline, y_vals_tangencyline, "k-");
        // }

        // 绘制凸多边形
        std::vector<double> rectangle_x, rectangle_y, polygon_x_simply, polygon_y_simply;
        std::string filename_rectangle = Safety_Corridor_Path + "Rectangle_Corridor_" + std::to_string(i) + ".csv";
        fm.read_csv(filename_rectangle, rectangle_x, rectangle_y);
        fm.read_csv(filename_vertices_no_iteration, polygon_x_simply, polygon_y_simply);
        rec_cor_area.push_back(polygonArea(rectangle_x, rectangle_y));
        cor_area_no_iteration.push_back(polygonArea(polygon_x_simply, polygon_y_simply));
        cor_area.push_back(polygonArea(poly_x, poly_y));

        for (int i = 0; i < obstacles.size() - 2; i++)
        {
            // 绘制障碍物
            if (i == 0)
            {
                plt::named_plot("Obstacle", obstacles[i].x, obstacles[i].y, "r-"); // 障碍物
            }
            else
            {
                plt::plot(obstacles[i].x, obstacles[i].y, "r-");
            }
        }
        plt::named_plot("Max Inscribed Ellipse", c_x_vals, c_y_vals, "k-");
        plt::named_plot("Convex Polygon Safe Corridor", poly_x, poly_y, "g-");
        plt::named_plot("Lattice Course", x_lattice, y_lattice, "b--");
        plt::named_plot("Safe Corridor with no optimization", polygon_x_simply, polygon_y_simply, "m--");
        plt::named_plot("Rectangle Safety Corridor", rectangle_x, rectangle_y, "y-");
        plt::named_plot("Vehicle", han.corner_points_x, han.corner_points_y, "k--"); // 黑色线条
        // 画椭圆中心
        plt::named_plot("center of Ellipse", std::vector<double>{center(0)}, std::vector<double>{center(1)}, "bo");
        plt::pause(0.2);
        plt::axis("equal");
        plt::legend();
        // 坐标轴
        plt::grid(true);
        plt::xlabel("X Axis");
        plt::ylabel("Y Axis");
        plt::title("Safety Corridor generation for Lattice planning");
    }
    std::vector<double> step_num_vector(rec_cor_area.size());
    std::iota(step_num_vector.begin(), step_num_vector.end(), 1);
    double rec_cor_area_mean = std::accumulate(rec_cor_area.begin(), rec_cor_area.end(), 0.0) / rec_cor_area.size();
    double cor_area_no_iter_mean = std::accumulate(cor_area_no_iteration.begin(), cor_area_no_iteration.end(), 0.0) / cor_area_no_iteration.size();
    double cor_area_mean = std::accumulate(cor_area.begin(), cor_area.end(), 0.0) / cor_area.size();

    std::string file_area_comparison = "../Data/Area_comparison/";
    std::string area_comparison_no_iter = file_area_comparison + "corridor_area_no_iteration.csv";
    std::string area_comparison = file_area_comparison + "corridor_area_with_stop_iteration_threshold_" + std::to_string(params.stop_iteration_criterion * 100) + ".csv";
    std::string area_comparison_rec = file_area_comparison + "rectangle_corridor_area.csv";
    save_to_csv_xy(area_comparison_no_iter, step_num_vector, cor_area_no_iteration);
    save_to_csv_xy(area_comparison, step_num_vector, cor_area);
    save_to_csv_xy(area_comparison_rec, step_num_vector, rec_cor_area);

    std::cout << "rec_cor_area: " << rec_cor_area_mean << std::endl;
    std::cout << "cor_area_no_iter: " << cor_area_no_iter_mean << std::endl;
    std::cout << "cor_area: " << cor_area_mean << std::endl;
    std::cout << "corridor increase magnitude compared to rectangle: " << (cor_area_mean - rec_cor_area_mean) * 100 / rec_cor_area_mean << "%" << std::endl;
    std::cout << "corridor increase magnitude compared to no optimization: " << (cor_area_mean - cor_area_no_iter_mean) * 100 / cor_area_no_iter_mean << "%" << std::endl;
    plt::figure();
    plt::named_plot("rectangle corridor area", step_num_vector, rec_cor_area, "r-");
    plt::named_plot("corridor area", step_num_vector, cor_area, "b-");
    plt::named_plot("corridor area no iteration", step_num_vector, cor_area_no_iteration, "g-");
    plt::legend();
    plt::xlabel("trajectory Index");
    plt::ylabel("area");
    plt::show();
    return 0;
}
