// #include <fmt/core.h>

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "cubic_spline.hpp"
#include "matplotlibcpp.h"
#include "quartic_polynomial.hpp"
#include "quintic_polynomial.hpp"
#include "utils.hpp"
#include "intersect_collision.h"

using std::vector;
using namespace Eigen;
namespace plt = matplotlibcpp;
typedef Eigen::Vector2d Vec2;

// 道路宽度常量
constexpr double ROAD_WIDTH = 3.75;
// 道路采样步长常量
constexpr double ROAD_SAMPLE_STEP = 1;
// 目标速度常量，单位转换为 m/s
constexpr double TARGET_SPEED = 10 / 3.6;
// 速度采样步长常量，单位转换为 m/s
constexpr double SPEED_SAMPLE_STEP = 1 / 3.6;

// 时间步长常量
constexpr double T_STEP = 0.15;
// 扭矩惩罚系数常量
constexpr double K_JERK = 0.1;
// 时间惩罚系数常量
constexpr double K_TIME = 1.0;
// 速度差惩罚系数常量
constexpr double K_V_DIFF = 1.0;
// 偏移惩罚系数常量
constexpr double K_OFFSET = 1.5;
// 碰撞惩罚系数常量
constexpr double K_COLLISION = 1000;

// 最大速度常量，单位转换为 m/s
constexpr double MAX_SPEED = 50.0 / 3.6;
// 最大加速度常量
constexpr double MAX_ACCEL = 8.0;
// 最大曲率常量
constexpr double MAX_CURVATURE = 1.27;

// 表示路径的类，包含时间、代价、横向和纵向的位移、速度、加速度、扭矩等
class Path
{
public:
    vector<double> t;
    double cost = 0.0;

    vector<double> l;
    vector<double> l_v;
    vector<double> l_a;
    vector<double> l_jerk;

    vector<double> s;
    vector<double> s_v;
    vector<double> s_a;
    vector<double> s_jerk;

    vector<double> x;
    vector<double> y;
    vector<double> yaw;
    vector<double> ds;
    vector<double> curv;

    Path() {}
    ~Path() {}

    // 将路径从 SL 坐标系转换为 XY 坐标系
    void SL_2_XY(CubicSpline2D &ref_path);
    // 计算路径的航向角和曲率
    void calc_yaw_curv(void);
    // 比较两个路径对象的代价
    bool operator<(const Path &other) const { return cost < other.cost; }
};

// 将路径从 SL 坐标系转换为 XY 坐标系
void Path::SL_2_XY(CubicSpline2D &ref_path)
{
    x.clear();
    y.clear();

    for (size_t i = 0; i < s.size(); ++i)
    {
        if (s[i] > ref_path.s.back())
        {
            break;
        }

        Vector2d xy_ref = ref_path.calc_position(s[i]);
        double yaw = ref_path.calc_yaw(s[i]);
        double x_ref = xy_ref[0] + l[i] * cos(yaw + M_PI_2);
        double y_ref = xy_ref[1] + l[i] * sin(yaw + M_PI_2);

        x.push_back(x_ref);
        y.push_back(y_ref);
    }
}

// 计算路径的航向角和曲率
void Path::calc_yaw_curv(void)
{
    yaw.clear();
    curv.clear();
    ds.clear();

    for (size_t i = 0; i + 1 < x.size(); ++i)
    {
        double dx = x[i + 1] - x[i];
        double dy = y[i + 1] - y[i];
        ds.push_back(hypot(dx, dy));
        yaw.push_back(atan2(dy, dx));
    }

    if (yaw.empty())
    {
        return;
    }
    yaw.push_back(yaw.back());
    ds.push_back(ds.back());

    for (size_t i = 0; i + 1 < yaw.size(); ++i)
    {
        curv.emplace_back((yaw[i + 1] - yaw[i]) / ds[i]);
    }
}

// 根据给定的 x, y 坐标生成参考路径
vector<vector<double>> get_reference_line(vector<double> cx, vector<double> cy,
                                          CubicSpline2D &spline)
{
    vector<double> x;
    vector<double> y;
    for (size_t idx = 0; idx < cx.size(); idx += 3)
    {
        x.push_back(cx[idx]);
        y.push_back(cy[idx]);
    }

    vector<vector<double>> traj = CubicSpline2D::calc_spline_course(x, y, 0.1);
    spline = CubicSpline2D(x, y);

    return traj;
}

void L_to_G_for_vehicle(Eigen::Matrix<double, 2, 4> &vehicle, const utils::VehicleConfig &vc, const double &x_point, const double &y_point, const double &theta, const bool &add_safety_margin)
{
    if (add_safety_margin)
    {
        vehicle << vc.RF + vc.SM, vc.RF + vc.SM, -vc.RB - vc.SM, -vc.RB - vc.SM, -vc.W / 2 - vc.SM, vc.W / 2 + vc.SM, vc.W / 2 + vc.SM, -vc.W / 2 - vc.SM;
    }
    else
    {
        vehicle << vc.RF, vc.RF, -vc.RB, -vc.RB, -vc.W / 2, vc.W / 2, vc.W / 2, -vc.W / 2;
    }
    Eigen::Matrix2d rot;
    rot << cos(theta), -sin(theta), sin(theta), cos(theta);
    vehicle = rot * vehicle;
    vehicle += Vec2(x_point, y_point).replicate(1, 4);
}

// 检查路径是否与障碍物发生碰撞
bool is_path_collision(const Path &path, const utils::VehicleConfig &vc,
                       const vector<Polygon_simple> &obs)
{
    vector<double> x;
    vector<double> y;
    vector<double> yaw;
    for (size_t i = 0; i < path.x.size(); i += 3)
    {
        x.push_back(path.x[i]);
        y.push_back(path.y[i]);
        yaw.push_back(path.yaw[i]);
    }

    for (size_t i = 0; i < x.size(); ++i)
    {
        for (size_t k = 0; k < obs.size(); k++)
        {
            Eigen::Matrix<double, 2, 4> vehicle;
            bool add_safety_margin = true;
            L_to_G_for_vehicle(vehicle, vc, x[i], y[i], yaw[i], add_safety_margin);
            Polygon_simple vehicle_global;
            vehicle_global.x = vector<double>{vehicle(0, 0), vehicle(0, 1), vehicle(0, 2), vehicle(0, 3), vehicle(0, 0)};
            vehicle_global.y = vector<double>{vehicle(1, 0), vehicle(1, 1), vehicle(1, 2), vehicle(1, 3), vehicle(1, 0)};
            if (EdgeIntersector::collisioncheck(vehicle_global, obs[k]))
            {
                return true;
            }
        }
    }

    return false;
}

// 验证路径的可行性，包括速度、加速度和曲率是否在允许范围内
bool verify_path(const Path &path, const utils::VehicleConfig &vc,
                 const vector<Polygon_simple> &obs)
{
    for (size_t i = 0; i < path.s_v.size(); ++i)
    {
        if (path.s_v[i] > MAX_SPEED || abs(path.s_a[i]) > MAX_ACCEL ||
            abs(path.curv[i]) > MAX_CURVATURE || is_path_collision(path, vc, obs))
        {
            return false;
        }
    }

    return true;
}

// 生成候选路径
vector<Path> sampling_paths(double l0, double l0_v, double l0_a, double s0, double s0_v,
                            double s0_a, CubicSpline2D &ref_path, const utils::VehicleConfig &vc,
                            const vector<Polygon_simple> &obs)
{
    vector<Path> paths;

    for (double s1_v = TARGET_SPEED * 0.2; s1_v < TARGET_SPEED * 1.2; s1_v += TARGET_SPEED * 0.2)
    {
        for (double t1 = 4.5; t1 < 5.5; t1 += 0.2)
        {
            Path path_pre;
            QuarticPolynomial path_lon(s0, s0_v, s0_a, s1_v, 0.0, t1);

            for (double t = 0.0; t < t1; t += T_STEP)
            {
                path_pre.t.push_back(t);
                path_pre.s.push_back(path_lon.calc_point(t));
                path_pre.s_v.push_back(path_lon.calc_first_derivative(t));
                path_pre.s_a.push_back(path_lon.calc_second_derivative(t));
                path_pre.s_jerk.push_back(path_lon.calc_third_derivative(t));
            }

            for (double l1 = -2 * ROAD_WIDTH; l1 < 2 * ROAD_WIDTH; l1 += ROAD_SAMPLE_STEP)
            {
                Path path = path_pre;
                QuinticPolynomial path_lat(l0, l0_v, l0_a, l1, 0.0, 0.0, t1);

                for (double t : path_pre.t)
                {
                    path.l.push_back(path_lat.calc_point(t));
                    path.l_v.push_back(path_lat.calc_first_derivative(t));
                    path.l_a.push_back(path_lat.calc_second_derivative(t));
                    path.l_jerk.push_back(path_lat.calc_third_derivative(t));
                }

                path.SL_2_XY(ref_path);
                path.calc_yaw_curv();
                if (path.yaw.empty())
                {
                    continue;
                }

                double l_jerk_sum = 0.0;
                double s_jerk_sum = 0.0;
                double v_diff = abs(TARGET_SPEED - path.s_v.back());
                for (size_t i = 0; i < path.l_jerk.size(); ++i)
                {
                    l_jerk_sum += abs(path.l_jerk[i]);
                    s_jerk_sum += abs(path.s_jerk[i]);
                }
                path.cost = K_JERK * (l_jerk_sum + s_jerk_sum) + K_V_DIFF * v_diff +
                            K_TIME * t1 * 2 + K_OFFSET * abs(path.l.back());
                paths.emplace_back(path);
            }
        }
    }

    return paths;
}

// 从候选路径中提取最优路径
Path extract_optimal_path(vector<Path> &paths, const utils::VehicleConfig &vc,
                          const vector<Polygon_simple> &obs)
{
    Path path;
    if (paths.empty())
    {
        return path;
    }

    std::sort(paths.begin(), paths.end());
    for (Path &p : paths)
    {
        if (verify_path(p, vc, obs))
        {
            path = p;
            break;
        }
    }

    return path;
}

// 用于巡航场景的路径规划
Path lattice_planner(double l0, double l0_v, double l0_a, double s0, double s0_v, double s0_a,
                     CubicSpline2D &ref_path, const utils::VehicleConfig &vc,
                     const vector<Polygon_simple> &obs)
{
    vector<Path> paths = sampling_paths(l0, l0_v, l0_a, s0, s0_v, s0_a, ref_path, vc, obs);
    Path path = extract_optimal_path(paths, vc, obs);

    return path;
}

void draw_vehicle(const utils::VehicleConfig &vc, const double &x_point, const double &y_point, const double &theta, const bool &is_ego)
{
    // Polygon_simple vehicle_global = transfer_polygon_to_global_coordinate(vehicle, x_point, y_point, theta);
    Eigen::Matrix<double, 2, 4> vehicle;
    bool add_safety_margin = false;
    L_to_G_for_vehicle(vehicle, vc, x_point, y_point, theta, add_safety_margin);
    Polygon_simple vehicle_global;
    vehicle_global.x = vector<double>{vehicle(0, 0), vehicle(0, 1), vehicle(0, 2), vehicle(0, 3), vehicle(0, 0)};
    vehicle_global.y = vector<double>{vehicle(1, 0), vehicle(1, 1), vehicle(1, 2), vehicle(1, 3), vehicle(1, 0)};
    if (is_ego)
    {
        plt::plot(vehicle_global.x, vehicle_global.y, "--k");
    }
    else
    {
        plt::plot(vehicle_global.x, vehicle_global.y, "r-");
    }
}

void save_to_csv(const std::string &filename,
                 const std::vector<double> &x,
                 const std::vector<double> &y,
                 const std::vector<double> &theta)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // 写入表头
    file << "x,y,theta\n";

    // 写入每一行数据
    for (size_t i = 0; i < x.size(); i += 3)
    {
        file << x[i] << "," << y[i] << "," << theta[i] << "\n";
    }

    file.close();
    std::cout << "数据已写入到: " << filename << std::endl;
}

// 模拟巡航场景
void cruise_case(const utils::VehicleConfig &vc)
{
    // CruiseRoadLine cruise_line;
    utils::file_manager fm;
    std::vector<double> ax, ay, ath;
    Polygon_simple vehicle_obs;
    std::vector<Polygon_simple> obs, obstacles;
    CubicSpline2D spline;
    // 基础路径
    std::string base_path = "../Data/Path_with_triangle_obstacles/";

    // 读取 Predicted_trajectory
    std::string predicted_trajectory_file = base_path + "reference_track.csv";
    fm.read_csv_with_theta(predicted_trajectory_file, ax, ay, ath);

    // 读取 obs 数据
    for (int i = 0; i <= 19; ++i)
    {
        std::string obs_file = base_path + "obs_" + std::to_string(i) + ".csv";
        Polygon_simple obstacle;
        fm.read_csv(obs_file, obstacle.x, obstacle.y);
        obstacles.push_back(obstacle);
        obs.push_back(obstacle);
    }

    bool add_safety_margin = false;
    std::vector<int> index_of_obs_vehicle = {10, 20};

    double lateral_offset = 1.5;
    for (int i = 0; i < index_of_obs_vehicle.size(); ++i)
    {
        lateral_offset = lateral_offset * -1.0;
        Eigen::Matrix<double, 2, 4> vehicle_ob;
        L_to_G_for_vehicle(vehicle_ob, vc, ax[index_of_obs_vehicle[i]], ay[index_of_obs_vehicle[i]] + lateral_offset, ath[index_of_obs_vehicle[i]], add_safety_margin);
        vehicle_obs.x = vector<double>{vehicle_ob(0, 0), vehicle_ob(0, 1), vehicle_ob(0, 2), vehicle_ob(0, 3), vehicle_ob(0, 0)};
        vehicle_obs.y = vector<double>{vehicle_ob(1, 0), vehicle_ob(1, 1), vehicle_ob(1, 2), vehicle_ob(1, 3), vehicle_ob(1, 0)};
        obs.push_back(vehicle_obs);
    }

    // 检查路径是否成功读取
    if (ax.empty() || ay.empty())
    {
        std::cerr << "Failed to read path points from file." << std::endl;
        // return -1;
        return;
    }

    // double delta_ax_last = ax[ax.size() - 1] - ax[ax.size() - 2];
    // double delta_ay_last = ay[ay.size() - 1] - ay[ay.size() - 2];
    // for (int i = 0; i < 2; ++i)
    // {
    //     ax.push_back(ax.back() + delta_ax_last);
    //     ay.push_back(ay.back() + delta_ay_last);
    // }

    vector<vector<double>> traj = get_reference_line(ax, ay, spline);

    double l0 = 0.0;          // 当前横向位置 [m]
    double l0_v = 0.0;        // 当前横向速度 [m/s]
    double l0_a = 0.0;        // 当前横向加速度 [m/s]
    double s0 = 0.0;          // 当前纵向位置
    double s0_v = 10.0 / 3.6; // 当前速度 [m/s]
    double s0_a = 0.0;
    vector<double> vehicle_x_history, vehicle_y_history, vehicle_yaw_history;
    utils::TicToc timer;
    std::vector<double> calculation_time;
    while (true)
    {
        timer.tic();
        Path path = lattice_planner(l0, l0_v, l0_a, s0, s0_v, s0_a, spline, vc, obs);
        calculation_time.push_back(timer.toc());

        if (path.x.empty())
        {
            // fmt::print("No feasible path found!!\n");
            std::cout << "No feasible path found!!" << std::endl;
            break;
        }

        l0 = path.l[1];
        l0_v = path.l_v[1];
        l0_a = path.l_a[1];
        s0 = path.s[1];
        s0_v = path.s_v[1];
        s0_a = path.s_a[1];

        if (hypot(path.x[1] - traj[0].back(), path.y[1] - traj[1].back()) <= 0.5)
        {
            // fmt::print("Goal\n");
            std::cout << "Goal" << std::endl;
            break;
        }

        double dy = (path.yaw[2] - path.yaw[1]) / path.ds[1];
        double steer = utils::pi_2_pi(atan(1.2 * vc.WB * dy));
        vehicle_x_history.emplace_back(path.x[1]);
        vehicle_y_history.emplace_back(path.y[1]);
        vehicle_yaw_history.emplace_back(path.yaw[1]);
        plt::cla();
        // plt::named_plot("Candidate trajectories", paths[0].x, paths[0].y, "-c");
        // for (size_t i = 1; i < paths.size(); i += (paths.size() / 10))
        // {
        //     plt::plot(paths[i].x, paths[i].y, "-c");
        // }
        plt::named_plot("Reference trajectory", ax, ay, "--b");

        // plt::named_plot("Optimal trajectory", path.x, path.y, "-g");
        for (int i = 0; i < obstacles.size(); i++)
        {
            if (i == 0)
            {
                plt::named_plot("Obstacles", obstacles[i].x, obstacles[i].y, "-r");
            }
            else
            {
                plt::plot(obstacles[i].x, obstacles[i].y, {{"color", "r"}});
            }
        }

        for (int i = 0; i < index_of_obs_vehicle.size(); i++)
        {
            if (i == 0)
            {
                draw_vehicle(vc, ax[index_of_obs_vehicle[i]], ay[index_of_obs_vehicle[i]] - 1.5, ath[index_of_obs_vehicle[i]], false);
            }
            else
            {
                draw_vehicle(vc, ax[index_of_obs_vehicle[i]], ay[index_of_obs_vehicle[i]] + 1.5, ath[index_of_obs_vehicle[i]], false);
            }
        }
        plt::named_plot("lattice planner trajectory", vehicle_x_history, vehicle_y_history, "-g");
        for (int k = 0; k < vehicle_x_history.size(); k++)
        {
            draw_vehicle(vc, vehicle_x_history[k], vehicle_y_history[k], vehicle_yaw_history[k], true);
        }
        // for (int k = 0; k < vehicle_x_history.size(); k = k + 3)
        // {
        //     draw_vehicle(vc, vehicle_x_history[k], vehicle_y_history[k], vehicle_yaw_history[k]);
        // }
        // plt::title("Lattice Planner in Cruising Scene V[km/h]:" +
        //            std::to_string(s0_v * 3.6).substr(0, 4));
        plt::axis("equal");
        plt::legend();
        plt::pause(0.0001);
    }
    // save_to_csv("../Data/Path_with_triangle_obstacles/lattice_path.csv", vehicle_x_history, vehicle_y_history, vehicle_yaw_history);
    double average_calculation_time = std::accumulate(calculation_time.begin(), calculation_time.end(), 0.0) / calculation_time.size();
    std::cout << "Average calculation time: " << average_calculation_time << " [ms]" << std::endl;
    plt::show();
}

int main()
{
    utils::VehicleConfig vc;
    cruise_case(vc);

    return 0;
}