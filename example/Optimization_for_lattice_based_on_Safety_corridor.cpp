#include <iostream>
#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <osqp/osqp.h>
#include <limits>
#include <vector>
#include <cmath>
#include "matplotlibcpp.h"
#include "utils.hpp"
using namespace Eigen;
typedef Eigen::Vector2d Vec2;
namespace plt = matplotlibcpp;
const double l = 2.92;

casadi::MX f(const casadi::MX &x, const casadi::MX &u)
{
    auto f1 = x(3) * casadi::MX::cos(x(2));     // v*cos(theta)
    auto f2 = x(3) * casadi::MX::sin(x(2));     // v*sin(theta)
    auto f3 = x(3) / l * casadi::MX::tan(x(4)); // v/l*tan(delta)
    auto f4 = u(0);                             // a
    auto f5 = u(1);                             // omega
    return vertcat(f1, f2, f3, f4, f5);
    // return vertcat(x(1), u - x(1));
}

casadi::MX loc_2_glo(const casadi::MX &x_car, const casadi::MX &y_car, const casadi::MX &theta_car, const double &x_loc, const double &y_loc)
{
    casadi::MX x_glo = x_car + x_loc * casadi::MX::cos(theta_car) - y_loc * casadi::MX::sin(theta_car);
    casadi::MX y_glo = y_car + x_loc * casadi::MX::sin(theta_car) + y_loc * casadi::MX::cos(theta_car);
    return vertcat(x_glo, y_glo);
}

class Obstacle
{
public:
    std::vector<double> x, y;
};

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

class Tangency_lines_single_SC
{
public:
    std::vector<double> a, b, c;
};

void Cal_curvature_standard(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &k)
{
    int N = x.size();
    k.resize(N, 0.0);
    for (int i = 1; i < N - 1; ++i)
    {
        double dx1 = x[i] - x[i - 1];
        double dy1 = y[i] - y[i - 1];
        double dx2 = x[i + 1] - x[i];
        double dy2 = y[i + 1] - y[i];

        double ds1 = std::sqrt(dx1 * dx1 + dy1 * dy1) + 1e-6;
        double ds2 = std::sqrt(dx2 * dx2 + dy2 * dy2) + 1e-6;

        double x_prime = (x[i + 1] - x[i - 1]) / (ds1 + ds2);
        double y_prime = (y[i + 1] - y[i - 1]) / (ds1 + ds2);
        double x_double_prime = (dx2 / ds2 - dx1 / ds1) / ((ds1 + ds2) / 2);
        double y_double_prime = (dy2 / ds2 - dy1 / ds1) / ((ds1 + ds2) / 2);

        double numerator = x_prime * y_double_prime - y_prime * x_double_prime;
        double denominator = std::pow(x_prime * x_prime + y_prime * y_prime, 1.5) + 1e-6;

        k[i] = numerator / denominator;
    }
    k[0] = k[1];
    k[N - 1] = k[N - 2];
}

double compute_mean_abs(const std::vector<double> &v)
{
    if (v.empty())
        return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0,
                           [](double sum, double val)
                           { return sum + std::abs(val); }) /
           v.size();
}

void save_to_csv(const std::string &filename,
                 const std::vector<double> &x,
                 const std::vector<double> &y)
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

int main()
{
    const double LF = 3.88;
    const double LR = 1.12;
    const double W = 1.92;
    const double WB = 2.92;
    const double DT = 0.45;

    utils::file_manager fm;
    std::vector<double> ax, ay, ath, x_lattice, y_lattice, theta_lattice;
    Obstacle obs_vehicle, obs_vehicle1, obs_lower_bound, obs_upper_bound;
    std::vector<Tangency_lines_single_SC> Tangency_lines_all_SCs;
    std::string base_path = "../Data/Path_with_triangle_obstacles/";
    std::string lattice_path_file = base_path + "lattice_path.csv";
    std::string reference_path_file = base_path + "reference_track.csv";
    std::string Safety_corridor_files = base_path + "Safety_corridor/";
    std::vector<Obstacle> obstacles, Safety_corridor_vertices;
    fm.read_csv_with_theta(lattice_path_file, x_lattice, y_lattice, theta_lattice);
    fm.read_csv_with_theta(reference_path_file, ax, ay, ath);

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
    // 读取 obs 数据
    for (int i = 0; i <= 19; ++i)
    {
        std::string obs_file = base_path + "obs_" + std::to_string(i) + ".csv";
        Obstacle obstacle;
        fm.read_csv(obs_file, obstacle.x, obstacle.y);
        obstacles.push_back(obstacle);
    }
    obstacles.push_back(obs_vehicle);
    obstacles.push_back(obs_vehicle1);

    for (int i = 0; i <= 42; ++i)
    {
        std::string SFV_file = Safety_corridor_files + "Polygon_vertices_" + std::to_string(i) + ".csv";
        Obstacle SFV;
        fm.read_csv(SFV_file, SFV.x, SFV.y);
        Safety_corridor_vertices.push_back(SFV);

        std::string TLS_file = Safety_corridor_files + "Tangency_lines_" + std::to_string(i) + ".csv";
        Tangency_lines_single_SC TLS;
        fm.read_csv_with_theta(TLS_file, TLS.a, TLS.b, TLS.c);
        Tangency_lines_all_SCs.push_back(TLS);
    }

    utils::TicToc timer;
    int N = x_lattice.size();           // number of nodes
    casadi::Opti opti = casadi::Opti(); // Optimization problem
    casadi::Slice all;
    // ---- Optimization variables ---------
    casadi::MX X = opti.variable(5, N); // state trajectory
    casadi::MX P = opti.parameter(10);  // track parameters
    casadi::MX U = opti.variable(2, N); // control

    // casadi::MX Time = opti.variable();

    auto pos_X = X(0, all);
    auto pos_Y = X(1, all);
    auto Theta = X(2, all);
    auto speed = X(3, all);
    auto delta = X(4, all);
    // auto DT = Time / N;

    auto init_X = P(0);
    auto init_Y = P(1);
    auto init_Theta = P(2);
    auto init_speed = P(3);
    auto init_delta = P(4);
    auto final_X = P(5);
    auto final_Y = P(6);
    auto final_Theta = P(7);
    auto final_speed = P(8);
    auto final_delta = P(9);
    auto a = U(0, all);
    auto omega = U(1, all);

    // ---- Cost function -----
    casadi::MX obj = 0;
    double w_reduce = 2.0;
    double w_speed = 1.0;
    double target_speed = 10.0;
    for (int k = 0; k < N; k++)
    {
        obj += w_reduce * (a(k) * a(k) + speed(k) * omega(k) * omega(k) * speed(k)) * DT;
        obj += w_speed * (speed(k) - target_speed) * (speed(k) - target_speed) / 100 * DT;
        // obj += (a(k) * a(k) + std::pow(vp.WB,2)/(casadi::MX::tan(delta(k)) * casadi::MX::tan(delta(k))))*DT;
    }

    // ---- dynamic constraints --------
    for (int k = 0; k < N; ++k)
    {
        casadi::MX k1 = f(X(all, k), U(all, k));
        casadi::MX x_next = X(all, k) + DT * k1;
        if (k < N - 1)
        {
            opti.subject_to(X(all, k + 1) == x_next);
        }
    }
    double w_terminal = 2000.0;
    obj += w_terminal * (casadi::MX::pow(pos_X(N - 1) - final_X, 2) +
                         casadi::MX::pow(pos_Y(N - 1) - final_Y, 2) +
                         casadi::MX::pow(Theta(N - 1) - final_Theta, 2));
    // ---- path constraints -----------
    opti.subject_to(0 <= speed <= 15);               // track speed limit
    opti.subject_to(-2 <= a <= 2);                   // control is limited
    opti.subject_to(-M_PI / 6 <= delta <= M_PI / 6); // steering angle is limited
    // opti.subject_to(-M_PI / 12 <= omega <= M_PI / 12); // steering rate is limited
    // ---- initial conditions --------
    opti.subject_to(pos_X(0) == init_X);     // longitudinal position
    opti.subject_to(pos_Y(0) == init_Y);     // lateral position
    opti.subject_to(Theta(0) == init_Theta); // initial heading angle
    opti.subject_to(speed(0) == init_speed); // initial speed
    opti.subject_to(delta(0) == init_delta); // initial steering angle
    // ---- final conditions -----------
    // opti.subject_to(pos_X(N - 1) == final_X);     // longitudinal position
    // opti.subject_to(pos_Y(N - 1) == final_Y);     // lateral position
    // opti.subject_to(Theta(N - 1) == final_Theta); // final heading angle
    // opti.subject_to(speed(N - 1) == final_speed); // final speed
    // opti.subject_to(delta(N - 1) == final_delta); // final steering angle

    // ---- Obstacle constraints ----------
    for (int k = 0; k < N; k++)
    {
        casadi::MX point1, point2, point3, point4;
        point1 = loc_2_glo(pos_X(k), pos_Y(k), Theta(k), LF, -W / 2.0);
        point2 = loc_2_glo(pos_X(k), pos_Y(k), Theta(k), LF, W / 2.0);
        point3 = loc_2_glo(pos_X(k), pos_Y(k), Theta(k), -LR, W / 2.0);
        point4 = loc_2_glo(pos_X(k), pos_Y(k), Theta(k), -LR, -W / 2.0);
        for (int i = 0; i < Tangency_lines_all_SCs[k].a.size(); i++)
        {
            double f_center = Tangency_lines_all_SCs[k].a[i] * x_lattice[k] + Tangency_lines_all_SCs[k].b[i] * y_lattice[k] + Tangency_lines_all_SCs[k].c[i];
            auto f_point1 = Tangency_lines_all_SCs[k].a[i] * point1(0) + Tangency_lines_all_SCs[k].b[i] * point1(1) + Tangency_lines_all_SCs[k].c[i];
            auto f_point2 = Tangency_lines_all_SCs[k].a[i] * point2(0) + Tangency_lines_all_SCs[k].b[i] * point2(1) + Tangency_lines_all_SCs[k].c[i];
            auto f_point3 = Tangency_lines_all_SCs[k].a[i] * point3(0) + Tangency_lines_all_SCs[k].b[i] * point3(1) + Tangency_lines_all_SCs[k].c[i];
            auto f_point4 = Tangency_lines_all_SCs[k].a[i] * point4(0) + Tangency_lines_all_SCs[k].b[i] * point4(1) + Tangency_lines_all_SCs[k].c[i];
            opti.subject_to(f_center * f_point1 > 0);
            opti.subject_to(f_center * f_point2 > 0);
            opti.subject_to(f_center * f_point3 > 0);
            opti.subject_to(f_center * f_point4 > 0);
        }
    }

    // ---- solver option ------
    casadi::Dict opts;
    // opts["ipopt.linear_solver"] = "ma57"; // 或 "ma27"/"pardiso"
    // opts["ipopt.tol"] = 1e-6;             // 宽松些能明显加速
    // opts["ipopt.acceptable_tol"] = 1e-4;
    // opts["ipopt.max_iter"] = 200; // 防止尾部迭代拉长
    // opts["ipopt.sb"] = "yes";
    // opts["print_time"] = 0;
    opts["ipopt.print_level"] = 0; // 输出的详细程度: 0为无输出, 3为最详细
    opts["ipopt.sb"] = "yes";
    opts["print_time"] = 0;
    //  将选项传递给求解器
    opti.solver("ipopt", opts);

    // ---- set initial and terminal states -----
    opti.set_value(init_X, x_lattice[0]);
    opti.set_value(init_Y, y_lattice[0]);
    opti.set_value(init_Theta, theta_lattice[0]);
    opti.set_value(init_speed, 0);
    opti.set_value(init_delta, 0);
    opti.set_value(final_X, x_lattice[N - 1]);
    opti.set_value(final_Y, y_lattice[N - 1]);
    opti.set_value(final_Theta, theta_lattice[N - 1]);
    opti.set_value(final_speed, 0);
    opti.set_value(final_delta, 0);

    // ---- set initial values for optimization  variables-----
    for (int k = 0; k < N; k++)
    {
        opti.set_initial(pos_X(k), x_lattice[k]);
        opti.set_initial(pos_Y(k), y_lattice[k]);
        opti.set_initial(Theta(k), theta_lattice[k]);
        opti.set_initial(speed(k), target_speed);
        opti.set_initial(delta(k), 0);
        opti.set_initial(a(k), 0);
        opti.set_initial(omega(k), 0);
    }

    double calculation_time = 0;
    // ---- solve NLP ------
    timer.tic();
    casadi::OptiSol sol = opti.solve(); // actual solve
    calculation_time += timer.toc();
    std::cout << "Calculation time: " << calculation_time << "ms" << std::endl;
    casadi::DM X_val = sol.value(X); // 变量 X 的值（位置和速度）
    casadi::DM U_val = sol.value(U); // 变量 U 的值（控制/油门）
    std::vector<double> opti_pos_x, opti_pos_y, opti_theta, opti_speed, opti_delta, opti_a, opti_omega, opti_time, opti_radius;
    // plt::cla();
    // // 1. 绘制障碍物
    // for (int i = 0; i < obstacles.size(); i++)
    // {
    //     // 绘制障碍物
    //     if (i == 0)
    //     {
    //         plt::named_plot("Obstacle", obstacles[i].x, obstacles[i].y, "r-"); // 障碍物
    //     }
    //     else
    //     {
    //         plt::plot(obstacles[i].x, obstacles[i].y, "r-");
    //     }
    // }

    for (int k = 0; k < N; k++)
    {
        plt::cla();
        // 1. 绘制障碍物
        for (int i = 0; i < obstacles.size(); i++)
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
        opti_time.push_back(k * DT);
        opti_pos_x.push_back((double)X_val(0, k));
        opti_pos_y.push_back((double)X_val(1, k));
        opti_theta.push_back((double)X_val(2, k));
        opti_speed.push_back((double)X_val(3, k));
        opti_delta.push_back((double)X_val(4, k));
        opti_a.push_back((double)U_val(0, k));
        opti_omega.push_back((double)U_val(1, k));
        opti_radius.push_back(WB / std::tan(opti_delta[k]));
        simple_vehicle_with_position controlled_vehicle(LF, LR, W, opti_pos_x[k], opti_pos_y[k], opti_theta[k]);
        std::vector<double> controlled_vehicle_x, controlled_vehicle_y;
        controlled_vehicle_x = controlled_vehicle.corner_points_x;
        controlled_vehicle_x.push_back(controlled_vehicle_x[0]);
        controlled_vehicle_y = controlled_vehicle.corner_points_y;
        controlled_vehicle_y.push_back(controlled_vehicle_y[0]);
        if (k == 0)
        {
            plt::named_plot("controlled vehicle", controlled_vehicle_x, controlled_vehicle_y, "k--"); // 车辆
        }
        else
        {
            plt::plot(controlled_vehicle_x, controlled_vehicle_y, "k--"); // 车辆
        }

        // 绘制障碍物
        if (k == 0)
        {
            plt::named_plot("Safety Corridor", Safety_corridor_vertices[k].x, Safety_corridor_vertices[k].y, "g-"); // 障碍物
        }
        else
        {
            plt::plot(Safety_corridor_vertices[k].x, Safety_corridor_vertices[k].y, "g-");
        }
        plt::named_plot("Optimized course", opti_pos_x, opti_pos_y, "k-");
        plt::named_plot("lattice path", x_lattice, y_lattice, "b--");
        plt::axis("equal");
        plt::legend();
        plt::grid(true);
        plt::xlabel("X Axis");
        plt::ylabel("Y Axis");
        // plt::pause(0.2);
    }

    // plt::title("Optimized outcome");

    std::vector<double> opti_curvature, lattice_curvature;
    Cal_curvature_standard(x_lattice, y_lattice, lattice_curvature);
    Cal_curvature_standard(opti_pos_x, opti_pos_y, opti_curvature);
    double lattice_curvature_mean = compute_mean_abs(lattice_curvature);
    double opti_curvature_mean = compute_mean_abs(opti_curvature);
    plt::figure();
    plt::named_plot("Optimized curvature", opti_time, opti_curvature, "r-");
    plt::named_plot("Lattice curvature", opti_time, lattice_curvature, "b--");
    plt::legend();
    plt::xlabel("Time");
    plt::ylabel("Curvature");
    // plt::title("Curvature");
    std::cout << "Lattice curvature magnitude mean: " << lattice_curvature_mean << std::endl;
    std::cout << "Optimized curvature magnitude mean: " << opti_curvature_mean << std::endl;
    std::cout << "Curvature improved magnitude: " << (lattice_curvature_mean - opti_curvature_mean) * 100.0 / lattice_curvature_mean << "%" << std::endl;
    plt::show();
    std::string Curvature_comparison = "../Data/Curvature_comparison";
    std::string Optimized_curvature = Curvature_comparison + "/Optimized_curvature.csv";
    std::string Lattice_curvature = Curvature_comparison + "/Lattice_curvature.csv";
    save_to_csv(Optimized_curvature, opti_time, opti_curvature);
    save_to_csv(Lattice_curvature, opti_time, lattice_curvature);

    return 0;
}