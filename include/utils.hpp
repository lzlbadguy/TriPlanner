#pragma once
#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <climits>
#include <numeric>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

enum class Gear
{
    GEAR_DRIVE,
    GEAR_REVERSE
};

namespace utils
{

    template <typename T>
    int sign(T num)
    {
        if (num < 0)
        {
            return -1;
        }

        return 1;
    }

    template <typename T>
    Eigen::Matrix3d transformation_matrix2d(T x, T y, T theta)
    {
        Eigen::Matrix3d trans;
        trans << cos(theta), -sin(theta), x, sin(theta), cos(theta), y, 0, 0, 1;

        return trans;
    }

    template <typename T>
    Eigen::Matrix2d rotation_matrix2d(T theta)
    {
        Eigen::Matrix2d rotation;
        rotation << cos(theta), -sin(theta), sin(theta), cos(theta);

        return rotation;
    }

    template <typename T>
    double pi_2_pi(T theta)
    {
        while (theta > M_PI)
        {
            theta -= 2.0 * M_PI;
        }
        while (theta < -M_PI)
        {
            theta += 2.0 * M_PI;
        }

        return theta;
    }

    template <typename T>
    T max(std::vector<T> vec)
    {
        int size = vec.size();
        assert(size > 0);

        T ret = vec[0];
        for (int idx = 1; idx < size; ++idx)
        {
            if (vec[idx] > ret)
            {
                ret = vec[idx];
            }
        }

        return ret;
    }

    template <typename T>
    T min(std::vector<T> vec)
    {
        int size = vec.size();
        assert(size > 0);

        T ret = vec[0];
        for (int idx = 1; idx < size; ++idx)
        {
            if (vec[idx] < ret)
            {
                ret = vec[idx];
            }
        }

        return ret;
    }

    template <typename T>
    std::vector<T> diff(const std::vector<T> &vec)
    {
        std::vector<T> ret;
        for (size_t idx = 1; idx < vec.size(); ++idx)
        {
            ret.push_back(vec[idx] - vec[idx - 1]);
        }

        return ret;
    }

    template <typename T>
    std::vector<T> cumsum(std::vector<T> vec)
    {
        std::vector<T> output;
        T tmp = 0;
        for (size_t idx = 0; idx < vec.size(); ++idx)
        {
            tmp += vec[idx];
            output.push_back(tmp);
        }

        return output;
    }

    template <typename T>
    int search_index(std::vector<T> nums, T target)
    {
        int left = 0, right = nums.size() - 1;
        while (left <= right)
        {
            int mid = (right - left) / 2 + left;
            int num = nums[mid];
            if (num == target)
            {
                return mid;
            }
            else if (num > target)
            {
                right = mid - 1;
            }
            else
            {
                left = mid + 1;
            }
        }

        return -1;
    }

    template <typename T>
    double variance(const std::vector<T> &data)
    {
        if (data.empty())
        {
            return 0.0;
        }

        double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();

        double variance = 0.0;
        for (const double &value : data)
        {
            variance += pow(value - mean, 2);
        }
        variance /= data.size();

        return variance;
    }

    class TicToc
    {
    public:
        TicToc(void) { tic(); }

        void tic(void) { start = std::chrono::system_clock::now(); }

        double toc(void)
        {
            end = std::chrono::system_clock::now();
            std::chrono::duration<double, std::milli> elapsed_seconds = end - start;
            return elapsed_seconds.count();
        }

    private:
        std::chrono::time_point<std::chrono::system_clock> start, end;
    };

    class file_manager
    {
    public:
        // 读取 CSV 文件的类函数
        static void read_csv(const std::string &file_name, std::vector<double> &ax, std::vector<double> &ay);
        static void read_csv_with_theta(const std::string &file_name, std::vector<double> &ax, std::vector<double> &ay, std::vector<double> &theta);
    };

    // file_manage 类成员函数的实现
    inline void file_manager::read_csv(const std::string &file_name, std::vector<double> &ax, std::vector<double> &ay)
    {
        std::ifstream file(file_name);
        if (!file.is_open())
        {
            std::cerr << "Unable to open file: " << file_name << std::endl;
            return;
        }

        std::string line;
        // 跳过第一行标题
        std::getline(file, line);

        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string x_str, y_str;
            if (std::getline(ss, x_str, ',') && std::getline(ss, y_str, ','))
            {
                ax.push_back(std::stod(x_str)); // 转为 double
                ay.push_back(std::stod(y_str));
            }
        }
        file.close();
    }

    inline void file_manager::read_csv_with_theta(const std::string &file_name, std::vector<double> &ax, std::vector<double> &ay, std::vector<double> &ath)
    {
        std::ifstream file(file_name);
        if (!file.is_open())
        {
            std::cerr << "Unable to open file: " << file_name << std::endl;
            return;
        }

        std::string line;
        // 跳过第一行标题
        std::getline(file, line);

        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string x_str, y_str, theta_str;
            if (std::getline(ss, x_str, ',') && std::getline(ss, y_str, ',')&& std::getline(ss, theta_str, ','))
            {
                ax.push_back(std::stod(x_str)); // 转为 double
                ay.push_back(std::stod(y_str));
                ath.push_back(std::stod(theta_str));
            }
        }
        file.close();
    }

    class VehicleConfig
    {
    public:
        double RF = 3.88; // [m] distance from rear to vehicle front end of vehicle
        double RB = 1.12; // [m] distance from rear to vehicle back end of vehicle
        double W = 1.92;  // [m] width of vehicle
        double WD = 0.7;  // [m] distance between left-right wheels
        double WB = 2.92; // [m] Wheel base
        double TR = 0.44; // [m] Tyre radius
        double TW = 0.7;  // [m] Tyre width
        double MAX_STEER = M_PI / 4.0;
        double MAX_ACCEL = 2;
        double MAX_SPEED = 30.0 / 3.6;
        double MIN_SPEED = -20 / 3.6;
        double SM = 0.2; // [m] Safety margin
    };
} // namespace utils

#endif // __UTILS_HPP
