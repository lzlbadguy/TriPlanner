import casadi as ca
import numpy as np
import math
import matplotlib.pyplot as plt
import csv
import os
from scipy.spatial.transform import Rotation as R

# Global constants
l = 2.92
LF = 3.88
LR = 1.12
W = 1.92
WB = 2.92
DT = 0.45

class Obstacle:
    def __init__(self):
        self.x = []
        self.y = []

class SimpleVehicleWithPosition:
    def __init__(self, LF_, LR_, W_, x_, y_, theta_):
        self.LF = LF_
        self.LR = LR_
        self.W = W_
        self.x = x_
        self.y = y_
        self.theta = theta_
        
        # Rotation matrix
        self.R = np.array([[math.cos(theta_), -math.sin(theta_)],
                          [math.sin(theta_), math.cos(theta_)]])
        
        # Corner points in local coordinates
        self.corner_points = [
            np.array([LF_, -W_ / 2]),
            np.array([LF_, W_ / 2]),
            np.array([-LR_, W_ / 2]),
            np.array([-LR_, -W_ / 2])
        ]
        
        # Transform to global coordinates
        self.corner_points_x = []
        self.corner_points_y = []
        for i in range(4):
            point_global = self.R @ self.corner_points[i] + np.array([x_, y_])
            self.corner_points_x.append(point_global[0])
            self.corner_points_y.append(point_global[1])

class TangencyLinesSingleSC:
    def __init__(self):
        self.a = []
        self.b = []
        self.c = []

def f(x, u):
    """Dynamics function"""
    f1 = x[3] * ca.cos(x[2])     # v*cos(theta)
    f2 = x[3] * ca.sin(x[2])     # v*sin(theta)
    f3 = x[3] / l * ca.tan(x[4]) # v/l*tan(delta)
    f4 = u[0]                    # a
    f5 = u[1]                    # omega
    return ca.vertcat(f1, f2, f3, f4, f5)

def loc_2_glo(x_car, y_car, theta_car, x_loc, y_loc):
    """Convert local coordinates to global coordinates"""
    x_glo = x_car + x_loc * ca.cos(theta_car) - y_loc * ca.sin(theta_car)
    y_glo = y_car + x_loc * ca.sin(theta_car) + y_loc * ca.cos(theta_car)
    return ca.vertcat(x_glo, y_glo)

def cal_curvature_standard(x, y):
    """Calculate curvature of a path"""
    N = len(x)
    k = np.zeros(N)
    
    for i in range(1, N - 1):
        dx1 = x[i] - x[i - 1]
        dy1 = y[i] - y[i - 1]
        dx2 = x[i + 1] - x[i]
        dy2 = y[i + 1] - y[i]
        
        ds1 = math.sqrt(dx1 * dx1 + dy1 * dy1) + 1e-6
        ds2 = math.sqrt(dx2 * dx2 + dy2 * dy2) + 1e-6
        
        x_prime = (x[i + 1] - x[i - 1]) / (ds1 + ds2)
        y_prime = (y[i + 1] - y[i - 1]) / (ds1 + ds2)
        x_double_prime = (dx2 / ds2 - dx1 / ds1) / ((ds1 + ds2) / 2)
        y_double_prime = (dy2 / ds2 - dy1 / ds1) / ((ds1 + ds2) / 2)
        
        numerator = x_prime * y_double_prime - y_prime * x_double_prime
        denominator = math.pow(x_prime * x_prime + y_prime * y_prime, 1.5) + 1e-6
        
        k[i] = numerator / denominator
    
    k[0] = k[1]
    k[N - 1] = k[N - 2]
    
    return k

def compute_mean_abs(v):
    """Compute mean of absolute values"""
    if len(v) == 0:
        return 0.0
    return sum(abs(val) for val in v) / len(v)

def save_to_csv(filename, x, y):
    """Save data to CSV file"""
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['x', 'y'])  # Header
        for i in range(len(x)):
            writer.writerow([x[i], y[i]])
    print(f"数据已写入到: {filename}")

def read_csv(filename, x_col=0, y_col=1):
    """Read CSV file and return x and y columns"""
    x, y = [], []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            x.append(float(row[x_col]))
            y.append(float(row[y_col]))
    return x, y

def read_csv_with_theta(filename, x_col=0, y_col=1, theta_col=2):
    """Read CSV file and return x, y, and theta columns"""
    x, y, theta = [], [], []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            x.append(float(row[x_col]))
            y.append(float(row[y_col]))
            theta.append(float(row[theta_col]))
    return x, y, theta

def main():
    # File paths
    base_path = "../Data/Path_with_triangle_obstacles/"
    lattice_path_file = base_path + "lattice_path.csv"
    reference_path_file = base_path + "reference_track.csv"
    safety_corridor_files = base_path + "Safety_corridor/"
    
    # Read lattice path and reference track
    x_lattice, y_lattice, theta_lattice = read_csv_with_theta(lattice_path_file)
    ax, ay, ath = read_csv_with_theta(reference_path_file)
    
    # Create obstacles
    obstacles = []
    safety_corridor_vertices = []
    tangency_lines_all_scs = []
    
    # Create vehicle obstacles
    han_obs = SimpleVehicleWithPosition(LF, LR, W, ax[10], ay[10] - 1.5, ath[10])
    han_obs1 = SimpleVehicleWithPosition(LF, LR, W, ax[20], ay[20] + 1.5, ath[10])
    
    obs_vehicle = Obstacle()
    obs_vehicle.x = han_obs.corner_points_x.copy()
    obs_vehicle.x.append(obs_vehicle.x[0])
    obs_vehicle.y = han_obs.corner_points_y.copy()
    obs_vehicle.y.append(obs_vehicle.y[0])
    
    obs_vehicle1 = Obstacle()
    obs_vehicle1.x = han_obs1.corner_points_x.copy()
    obs_vehicle1.x.append(obs_vehicle1.x[0])
    obs_vehicle1.y = han_obs1.corner_points_y.copy()
    obs_vehicle1.y.append(obs_vehicle1.y[0])
    
    obstacles.append(obs_vehicle)
    obstacles.append(obs_vehicle1)
    
    # Read obstacle data
    for i in range(20):  # 0 to 19
        obs_file = base_path + f"obs_{i}.csv"
        if os.path.exists(obs_file):
            obstacle = Obstacle()
            obstacle.x, obstacle.y = read_csv(obs_file)
            obstacles.append(obstacle)
    
    # Read safety corridor data
    for i in range(43):  # 0 to 42
        # Read polygon vertices
        sfv_file = safety_corridor_files + f"Polygon_vertices_{i}.csv"
        if os.path.exists(sfv_file):
            sfv = Obstacle()
            sfv.x, sfv.y = read_csv(sfv_file)
            safety_corridor_vertices.append(sfv)
        
        # Read tangency lines
        tls_file = safety_corridor_files + f"Tangency_lines_{i}.csv"
        if os.path.exists(tls_file):
            tls = TangencyLinesSingleSC()
            tls.a, tls.b, tls.c = read_csv_with_theta(tls_file)
            tangency_lines_all_scs.append(tls)
    
    # Optimization problem setup
    N = len(x_lattice)  # number of nodes
    opti = ca.Opti()  # Optimization problem
    
    # Optimization variables
    X = opti.variable(5, N)  # state trajectory
    P = opti.parameter(10)   # track parameters
    U = opti.variable(2, N)  # control
    
    # State variables
    pos_X = X[0, :]
    pos_Y = X[1, :]
    Theta = X[2, :]
    speed = X[3, :]
    delta = X[4, :]
    
    # Parameter indices
    init_X = P[0]
    init_Y = P[1]
    init_Theta = P[2]
    init_speed = P[3]
    init_delta = P[4]
    final_X = P[5]
    final_Y = P[6]
    final_Theta = P[7]
    final_speed = P[8]
    final_delta = P[9]
    
    # Control variables
    a = U[0, :]
    omega = U[1, :]
    
    # Cost function
    obj = 0
    w_reduce = 2.0
    w_speed = 1.0
    target_speed = 10.0
    
    for k in range(N):
        obj += w_reduce * (a[k]**2 + speed[k] * omega[k]**2 * speed[k]) * DT
        obj += w_speed * (speed[k] - target_speed)**2 / 100 * DT
    
    # Terminal cost
    w_terminal = 2000.0
    obj += w_terminal * ((pos_X[N-1] - final_X)**2 +
                         (pos_Y[N-1] - final_Y)**2 +
                         (Theta[N-1] - final_Theta)**2)
    
    # Dynamic constraints
    for k in range(N):
        k1 = f(X[:, k], U[:, k])
        x_next = X[:, k] + DT * k1
        if k < N - 1:
            opti.subject_to(X[:, k+1] == x_next)
    
    # Path constraints
    opti.subject_to(opti.bounded(0, speed, 15))              # track speed limit
    opti.subject_to(opti.bounded(-2, a, 2))                 # control is limited
    opti.subject_to(opti.bounded(-math.pi/6, delta, math.pi/6))  # steering angle is limited
    
    # Initial conditions
    opti.subject_to(pos_X[0] == init_X)      # longitudinal position
    opti.subject_to(pos_Y[0] == init_Y)      # lateral position
    opti.subject_to(Theta[0] == init_Theta)  # initial heading angle
    opti.subject_to(speed[0] == init_speed)  # initial speed
    opti.subject_to(delta[0] == init_delta)  # initial steering angle
    
    # Obstacle constraints
    for k in range(N):
        point1 = loc_2_glo(pos_X[k], pos_Y[k], Theta[k], LF, -W / 2.0)
        point2 = loc_2_glo(pos_X[k], pos_Y[k], Theta[k], LF, W / 2.0)
        point3 = loc_2_glo(pos_X[k], pos_Y[k], Theta[k], -LR, W / 2.0)
        point4 = loc_2_glo(pos_X[k], pos_Y[k], Theta[k], -LR, -W / 2.0)
        
        for i in range(len(tangency_lines_all_scs[k].a)):
            f_center = tangency_lines_all_scs[k].a[i] * x_lattice[k] + \
                        tangency_lines_all_scs[k].b[i] * y_lattice[k] + \
                        tangency_lines_all_scs[k].c[i]
            f_point1 = tangency_lines_all_scs[k].a[i] * point1[0] + \
                         tangency_lines_all_scs[k].b[i] * point1[1] + \
                         tangency_lines_all_scs[k].c[i]
            f_point2 = tangency_lines_all_scs[k].a[i] * point2[0] + \
                         tangency_lines_all_scs[k].b[i] * point2[1] + \
                         tangency_lines_all_scs[k].c[i]
            f_point3 = tangency_lines_all_scs[k].a[i] * point3[0] + \
                         tangency_lines_all_scs[k].b[i] * point3[1] + \
                         tangency_lines_all_scs[k].c[i]
            f_point4 = tangency_lines_all_scs[k].a[i] * point4[0] + \
                         tangency_lines_all_scs[k].b[i] * point4[1] + \
                         tangency_lines_all_scs[k].c[i]
            opti.subject_to(f_center * f_point1 > 0)
            opti.subject_to(f_center * f_point2 > 0)
            opti.subject_to(f_center * f_point3 > 0)
            opti.subject_to(f_center * f_point4 > 0)
    
    # Solver options
    opts = {}
    opts["ipopt.print_level"] = 0  # Output level: 0 for no output, 3 for most detailed
    opts["ipopt.sb"] = "yes"
    opts["print_time"] = 0
    opti.solver("ipopt", opts)
    
    # Set initial and terminal states
    opti.set_value(init_X, x_lattice[0])
    opti.set_value(init_Y, y_lattice[0])
    opti.set_value(init_Theta, theta_lattice[0])
    opti.set_value(init_speed, 0)
    opti.set_value(init_delta, 0)
    opti.set_value(final_X, x_lattice[N-1])
    opti.set_value(final_Y, y_lattice[N-1])
    opti.set_value(final_Theta, theta_lattice[N-1])
    opti.set_value(final_speed, 0)
    opti.set_value(final_delta, 0)
    
    # Set initial values for optimization variables
    for k in range(N):
        opti.set_initial(pos_X[k], x_lattice[k])
        opti.set_initial(pos_Y[k], y_lattice[k])
        opti.set_initial(Theta[k], theta_lattice[k])
        opti.set_initial(speed[k], target_speed)
        opti.set_initial(delta[k], 0)
        opti.set_initial(a[k], 0)
        opti.set_initial(omega[k], 0)
    
    # Solve NLP
    sol = opti.solve()  # actual solve
    
    X_val = sol.value(X)  # Variable X values (position and speed)
    U_val = sol.value(U)  # Variable U values (control/throttle)
    
    # Extract solution
    opti_pos_x, opti_pos_y, opti_theta = [], [], []
    opti_speed, opti_delta, opti_a, opti_omega = [], [], [], []
    opti_time, opti_radius = [], []
    
    for k in range(N):
        opti_time.append(k * DT)
        opti_pos_x.append(float(X_val[0, k]))
        opti_pos_y.append(float(X_val[1, k]))
        opti_theta.append(float(X_val[2, k]))
        opti_speed.append(float(X_val[3, k]))
        opti_delta.append(float(X_val[4, k]))
        opti_a.append(float(U_val[0, k]))
        opti_omega.append(float(U_val[1, k]))
        opti_radius.append(WB / math.tan(opti_delta[k]) if math.tan(opti_delta[k]) != 0 else float('inf'))
    
    # Calculate curvature
    opti_curvature = cal_curvature_standard(opti_pos_x, opti_pos_y)
    lattice_curvature = cal_curvature_standard(x_lattice, y_lattice)
    
    # Compute mean curvature
    lattice_curvature_mean = compute_mean_abs(lattice_curvature)
    opti_curvature_mean = compute_mean_abs(opti_curvature)
    
    # Print results
    print(f"Lattice curvature magnitude mean: {lattice_curvature_mean}")
    print(f"Optimized curvature magnitude mean: {opti_curvature_mean}")
    improvement = (lattice_curvature_mean - opti_curvature_mean) * 100.0 / lattice_curvature_mean
    print(f"Curvature improved magnitude: {improvement}%")
    
    # Save to CSV
    curvature_comparison = "../Data/Path_with_triangle_obstacles/Curvature_comparison"
    optimized_curvature = curvature_comparison + "/Optimized_curvature.csv"
    lattice_curvature_file = curvature_comparison + "/Lattice_curvature.csv"
    save_to_csv(optimized_curvature, opti_time, opti_curvature)
    save_to_csv(lattice_curvature_file, opti_time, lattice_curvature)
    
    # Visualization
    plt.figure(figsize=(10, 8))
    
    # Plot obstacles
    for i, obstacle in enumerate(obstacles):
        if i == 0:
            plt.plot(obstacle.x, obstacle.y, 'r-', label='Obstacle')
        else:
            plt.plot(obstacle.x, obstacle.y, 'r-')
    
    # Plot safety corridors and optimized path
    for k in range(N):
        # Plot safety corridor
        # if k < len(safety_corridor_vertices):
        #     if k == 0:
        #         plt.plot(safety_corridor_vertices[k].x, safety_corridor_vertices[k].y, 'g-', label='Safety Corridor')
        #     else:
        #         plt.plot(safety_corridor_vertices[k].x, safety_corridor_vertices[k].y, 'g-')
        
        # Plot vehicle positions
        if k % 1 == 0:  # Plot every 3rd vehicle position
            controlled_vehicle = SimpleVehicleWithPosition(
                LF, LR, W, opti_pos_x[k], opti_pos_y[k], opti_theta[k])
            controlled_vehicle_x = controlled_vehicle.corner_points_x + [controlled_vehicle.corner_points_x[0]]
            controlled_vehicle_y = controlled_vehicle.corner_points_y + [controlled_vehicle.corner_points_y[0]]
            plt.plot(controlled_vehicle_x, controlled_vehicle_y, 'k--', alpha=0.7)
    
    # Plot paths
    plt.plot(opti_pos_x, opti_pos_y, 'k-', label='Optimized course')
    plt.plot(x_lattice, y_lattice, 'b--', label='Lattice path')
    
    plt.axis('equal')
    plt.legend()
    plt.grid(True)
    plt.xlabel('X Axis')
    plt.ylabel('Y Axis')
    plt.title('Optimized outcome')
    plt.show()
    
    # Plot curvature
    plt.figure()
    plt.plot(opti_time, opti_curvature, 'r-', label='Optimized curvature')
    plt.plot(opti_time, lattice_curvature, 'b--', label='Lattice curvature')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Curvature')
    plt.title('Curvature')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()