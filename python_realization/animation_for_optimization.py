import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import csv
import math

# 创建目录用于保存动画
output_dir = "../Data/Path_with_triangle_obstacles/Animation"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 定义类
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

# 文件路径
base_path = "../Data/Path_with_triangle_obstacles/"
lattice_path_file = base_path + "lattice_path.csv"
reference_path_file = base_path + "reference_track.csv"
optimized_path_file = base_path + "optimized_path.csv"  # 优化后的路径文件
safety_corridor_files = base_path + "Safety_corridor/"

# 读取路径数据
x_lattice, y_lattice, theta_lattice = read_csv_with_theta(lattice_path_file)
ref_x, ref_y, ref_theta = read_csv_with_theta(reference_path_file)

# 读取优化后的路径数据
opti_pos_x, opti_pos_y, opti_theta = read_csv_with_theta(optimized_path_file)

# 创建障碍物
obstacles = []
safety_corridor_vertices = []

# 创建车辆障碍物
LF, LR, W = 3.88, 1.12, 1.92
han_obs = SimpleVehicleWithPosition(LF, LR, W, ref_x[10], ref_y[10] - 1.5, ref_theta[10])
han_obs1 = SimpleVehicleWithPosition(LF, LR, W, ref_x[20], ref_y[20] + 1.5, ref_theta[10])

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

# 读取障碍物数据
for i in range(20):  # 0 to 19
    obs_file = base_path + f"obs_{i}.csv"
    if os.path.exists(obs_file):
        obstacle = Obstacle()
        obstacle.x, obstacle.y = read_csv(obs_file)
        obstacles.append(obstacle)

# 读取安全走廊数据
for i in range(min(43, len(x_lattice))):  # 读取所有43个安全走廊
    sfv_file = safety_corridor_files + f"Polygon_vertices_{i}.csv"
    if os.path.exists(sfv_file):
        sfv = Obstacle()
        sfv.x, sfv.y = read_csv(sfv_file)
        safety_corridor_vertices.append(sfv)

# 创建图形和轴
fig, ax_plot = plt.subplots(figsize=(12, 10))
ax_plot.set_xlim(min(min(opti_pos_x), min(x_lattice)) - 5, max(max(opti_pos_x), max(x_lattice)) + 5)
ax_plot.set_ylim(min(min(opti_pos_y), min(y_lattice)) - 5, max(max(opti_pos_y), max(y_lattice)) + 5)
ax_plot.set_xlabel('X Axis')
ax_plot.set_ylabel('Y Axis')
ax_plot.set_title('Optimization Process Animation')
ax_plot.axis('equal')
ax_plot.grid(True)

# 绘制静态元素
# 绘制障碍物
for i, obstacle in enumerate(obstacles):
    if i == 0:
        ax_plot.plot(obstacle.x, obstacle.y, 'r-', label='Obstacle', linewidth=2)
    else:
        ax_plot.plot(obstacle.x, obstacle.y, 'r-', linewidth=2)

# 绘制晶格路径
ax_plot.plot(x_lattice, y_lattice, 'b--', label='Lattice path', linewidth=1.5)

# 绘制参考轨迹
# ax_plot.plot(ref_x, ref_y, 'g-', label='Reference track', linewidth=1.5)

# 初始化动态元素
path_line, = ax_plot.plot([], [], 'k-', linewidth=2, label='Optimized path')
vehicle_shape, = ax_plot.plot([], [], 'k--', linewidth=1.5, alpha=0.7, label='Vehicle position')
current_point, = ax_plot.plot([], [], 'ko', markersize=6, label='Current position')
safety_corridor_shape, = ax_plot.plot([], [], 'g-', linewidth=2, alpha=0.7, label='Safety corridor')

# 添加图例
ax_plot.legend(loc='upper right')

# 初始化函数：绘制背景
def init():
    path_line.set_data([], [])
    vehicle_shape.set_data([], [])
    current_point.set_data([], [])
    safety_corridor_shape.set_data([], [])
    return path_line, vehicle_shape, current_point, safety_corridor_shape

# 动画函数：在每一帧中更新图形
def animate(frame):
    # 更新优化路径
    path_line.set_data(opti_pos_x[:frame+1], opti_pos_y[:frame+1])
    
    # 更新当前位置点
    if frame < len(opti_pos_x):
        current_point.set_data([opti_pos_x[frame]], [opti_pos_y[frame]])
        
        # 更新车辆形状（根据opti_theta[frame]计算）
        vehicle = SimpleVehicleWithPosition(LF, LR, W, opti_pos_x[frame], opti_pos_y[frame], opti_theta[frame])
        vehicle_x = vehicle.corner_points_x + [vehicle.corner_points_x[0]]
        vehicle_y = vehicle.corner_points_y + [vehicle.corner_points_y[0]]
        vehicle_shape.set_data(vehicle_x, vehicle_y)
        
        # 更新安全走廊形状
        if frame < len(safety_corridor_vertices):
            corridor = safety_corridor_vertices[frame]
            if corridor.x and corridor.y:  # 确保数据存在
                corridor_x = corridor.x + [corridor.x[0]] if corridor.x[0] != corridor.x[-1] else corridor.x
                corridor_y = corridor.y + [corridor.y[0]] if corridor.y[0] != corridor.y[-1] else corridor.y
                safety_corridor_shape.set_data(corridor_x, corridor_y)
            else:
                safety_corridor_shape.set_data([], [])
        else:
            safety_corridor_shape.set_data([], [])
    
    return path_line, vehicle_shape, current_point, safety_corridor_shape

# 创建动画对象
N = len(opti_pos_x)
ani = animation.FuncAnimation(fig, animate, frames=N, init_func=init, blit=True, interval=200, repeat=True)

# 保存为GIF文件
print("正在生成动画...")
try:
    ani.save(os.path.join(output_dir, 'optimization_animation.gif'), writer='pillow', fps=5, dpi=100)
    print(f"动画已保存到: {os.path.join(output_dir, 'optimization_animation.gif')}")
except Exception as e:
    print(f"保存动画时出错: {e}")
    print("无法保存动画文件")

# 显示动画
plt.tight_layout()
plt.show()

print("动画生成完成！")

