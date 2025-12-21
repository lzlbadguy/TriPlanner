<div align="center">

# TriPlanner

[English](./README.md) | [中文](./README.zh-CN.md)

</div>

This repository implements the core ideas of the paper  
**“TriPlanner: A Tri-stage Planner for Safe Vehicle Navigation in Complex Environments”**  
(see the PDF in the root directory).

TriPlanner adopts a **three-stage planning pipeline** to generate safe and feasible vehicle trajectories in complex environments:
1. lattice-based candidate trajectory generation,
2. real-time safety corridor construction,
3. trajectory optimization within the corridor.

The project is implemented in **C++20** and integrates **CasADi**, **Eigen**, and related libraries for numerical computation, constrained optimization, and visualization.

---

## Three-Stage Planning Pipeline

- **Stage 1: Lattice Trajectory Generation**  
  Candidate trajectories are generated in the Frenet frame using quartic and quintic polynomials.  
  Constraints on velocity, curvature, and collision are checked, and feasible candidates are transformed into the global frame to serve as initial trajectories for subsequent stages.

- **Stage 2: Safety Corridor Generation**  
  Based on obstacle geometry and candidate trajectories, piecewise feasible safety corridors are constructed using **minimum-volume inscribed ellipses (MVIE)** and tangency constraints.  
  The resulting convex polygons provide explicit feasible regions for optimization.

- **Stage 3: Trajectory Optimization within the Corridor**  
  A nonlinear optimization problem is formulated under vehicle dynamics and corridor constraints using **CasADi + Ipopt**, producing smooth, dynamically feasible, and safe trajectories.  
  Visualization and data export are also supported.

---

## Repository Structure

- `example/`: Three example programs corresponding to the three planning stages in the paper.
  - `lattice_planner_optimization_pre_generation.cpp`  
    Stage 1: lattice candidate trajectory generation and filtering.
  - `Real_Time_Corridor_generation_for_lattice_based_Trajectory.cpp`  
    Stage 2: real-time safety corridor construction based on obstacles and trajectories.
  - `Optimization_for_lattice_based_on_Safety_corridor.cpp`  
    Stage 3: trajectory optimization with vehicle dynamics and safety corridor constraints.

- `include/`: Geometric and optimization components, including spline interpolation, minimum-volume ellipse computation, safety region solvers, collision checking, and a wrapper for `matplotlibcpp`.

- `Data/`: Reference trajectories, obstacle configurations, and safety corridor data used in the experiments.

- `CMakeLists.txt`: Build configuration for generating the example executables.

---

## Dependencies

- **Build**:  
  CMake ≥ 3.12, a C++20 compiler (gcc/g++ or clang/clang++), Threads.

- **Math & Optimization**:  
  Eigen3, CasADi (typically deployed with Ipopt).  
  Implementations related to minimum-volume ellipses and safety corridors are located in `include/`.

- **Visualization**:  
  `matplotlibcpp` (provided in `include/`), which depends on Python3 development headers and the Python `matplotlib` package.

- **Others**:  
  Python3 development headers (`Python3::Python`), standard C++ STL.

---

## CasADi and Ipopt Installation

Please refer to:  
[CasADi&Ipopt Install Instructions](./Casadi_Ipopt_install_instruction.md)
---

## Build and Run

```bash
mkdir -p build && cd build
cmake ..
make -j
## Executables are generated in build/, named after the example source files
./lattice_planner_optimization_pre_generation
./Real_Time_Corridor_generation_for_lattice_based_Trajectory
./Optimization_for_lattice_based_on_Safety_corridor
```