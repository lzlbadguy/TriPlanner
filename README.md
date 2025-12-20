# TriPlanner

本仓库实现了论文《TriPlanner: A Tri-stage Planner for Safe Vehicle Navigation in Complex Environments》（见根目录 PDF）的核心思路：以“三阶段”方式完成车辆在复杂场景中的安全轨迹规划——先生成格点候选，再实时构造安全通道，最后在通道内进行优化。项目以 C++20 编写，并结合 CasADi、Eigen 等库进行数值计算、约束优化与可视化。

## 三阶段规划流程
- **阶段 1：格点轨迹生成**（lattice generation）  
  在 Frenet 坐标系中用四次/五次多项式生成大量候选路径，检查速度、曲率、碰撞等约束，并转换到全局坐标，得到供后续使用的初始轨迹。
- **阶段 2：安全通道构造**（safety corridor generation）  
  结合障碍物轮廓与候选轨迹，使用最小体积内接椭圆（MVIE）和切线约束构造分段可行的安全多边形/通道，确保后续优化有明确可行域。
- **阶段 3：通道内优化**（trajectory optimization within corridor）  
  在车辆动力学模型与安全通道约束下搭建非线性优化（CasADi&Ipopt），得到满足动力学、安全和舒适度的平滑轨迹，并支持可视化和导出。

## 仓库结构
- `example/`: 三个示例程序，对应论文的三阶段流程。
  - `lattice_planner_optimization_pre_generation.cpp`: 阶段 1，生成与筛选格点候选轨迹。
  - `Real_Time_Corridor_generation_for_lattice_based_Trajectory.cpp`: 阶段 2，根据障碍物和轨迹实时生成安全通道。
  - `Optimization_for_lattice_based_on_Safety_corridor.cpp`: 阶段 3，将安全通道与车辆动力学输入优化器，输出平滑可行轨迹。
- `include/`: 规划所需的几何与优化组件（样条插值、最小体积椭圆、安全区域求解、碰撞检测、matplotlibcpp 封装等）。
- `Data/`: 论文实验使用的参考轨迹、障碍物、安全通道数据示例。
- `CMakeLists.txt`: 构建配置，生成各示例可执行文件。

## 依赖项
- 构建：CMake ≥ 3.12，C++20 编译器（gcc/g++ 或 clang/clang++），Threads。
- 数学/优化：Eigen3，CasADi（通常随 Ipopt 部署）,最小体积椭圆与安全通道相关的实现位于 `include/`。
- 可视化：`matplotlibcpp`（随 `include` 提供）依赖 Python3 头文件及 `matplotlib` Python 包。
- 其他：Python3 开发头文件（`Python3::Python`），标准 C++ STL。

## CasADi 与 Ipopt 部署：请参考 https://blog.csdn.net/qq_45090497/article/details/138240039

## 构建与运行
```bash
mkdir -p build && cd build
cmake ..
make -j
# 可执行文件位于 build/，与 example 源文件同名
```

运行示例（确保工作目录能访问 `Data/` 相对路径）：
```bash
./lattice_planner_optimization_pre_generation
./Real_Time_Corridor_generation_for_lattice_based_Trajectory
./Optimization_for_lattice_based_on_Safety_corridor
```
