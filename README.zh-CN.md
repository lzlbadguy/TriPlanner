


# TriPlanner

[English](./README.md) | [中文](./README.zh-CN.md)

</div>

本仓库实现了论文  
**《TriPlanner: A Tri-stage Planner for Safe Vehicle Navigation in Complex Environments》**  
（见根目录 PDF）的核心思路。

TriPlanner 采用**三阶段规划框架**，在复杂环境中完成车辆安全轨迹规划：
1. 格点候选轨迹生成，
2. 实时安全通道构造，
3. 通道内轨迹优化。

项目使用 **C++20** 编写，并结合 **CasADi、Eigen** 等库进行数值计算、约束优化与可视化。

---

## 三阶段规划流程

- **阶段 1：格点轨迹生成（Lattice Generation）**  
  在 Frenet 坐标系中使用四次/五次多项式生成大量候选轨迹，对速度、曲率和碰撞约束进行筛选，并转换到全局坐标系，得到后续阶段使用的初始轨迹。

- **阶段 2：安全通道构造（Safety Corridor Generation）**  
  结合障碍物轮廓与候选轨迹，基于**最大体积内接椭圆（MVIE）**和切线约束构造分段安全多边形通道，为优化阶段提供明确的可行域。

- **阶段 3：通道内优化（Trajectory Optimization within Corridor）**  
  在车辆动力学模型和安全通道约束下，使用 **CasADi + Ipopt** 构建非线性优化问题，生成平滑、可行且安全的最终轨迹，并支持可视化与数据导出。

---

## 仓库结构

- `example/`：三个示例程序，对应论文中的三阶段流程。
  - `lattice_planner_optimization_pre_generation.cpp`：阶段 1，格点候选轨迹生成与筛选。
  - `Real_Time_Corridor_generation_for_lattice_based_Trajectory.cpp`：阶段 2，实时安全通道生成。
  - `Optimization_for_lattice_based_on_Safety_corridor.cpp`：阶段 3，基于安全通道的轨迹优化。

- `include/`：规划所需的几何与优化组件，包括样条插值、最大体积椭圆计算、安全区域求解、碰撞检测以及 `matplotlibcpp` 封装。

- `Data/`：论文实验中使用的参考轨迹、障碍物与安全通道示例数据。

- `CMakeLists.txt`：构建配置文件，用于生成各示例可执行程序。

---

## 依赖项

- **构建**：  
  CMake ≥ 3.12，支持 C++20 的编译器（gcc/g++ 或 clang/clang++），Threads。

- **数学与优化**：  
  Eigen3，CasADi（通常与 Ipopt 一同部署）。  
  与最大体积椭圆和安全通道相关的实现位于 `include/` 目录。

- **可视化**：  
  `matplotlibcpp`（随 `include` 提供），依赖 Python3 开发头文件及 Python `matplotlib` 包。

- **其他**：  
  Python3 开发头文件（`Python3::Python`），标准 C++ STL。

---

## CasADi 与 Ipopt 部署

请参考：  
https://blog.csdn.net/qq_45090497/article/details/138240039

---

## 构建与运行

```bash
mkdir -p build && cd build
cmake ..
make -j
# 可执行文件位于 build/，名称与 example 中源文件一致
./lattice_planner_optimization_pre_generation
./Real_Time_Corridor_generation_for_lattice_based_Trajectory
./Optimization_for_lattice_based_on_Safety_corridor
```