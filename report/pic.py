import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

# 创建一个图形和坐标轴
fig, ax = plt.subplots()

# 设置绘图范围
ax.set_xlim(-0.1, 2.1)
ax.set_ylim(-0.1, 2.1)

# 网格参数
nx, ny = 3, 2  # 网格点的数量
dx, dy = 1, 1  # 网格间距

# 生成网格点
x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, 3)
X, Y = np.meshgrid(x, y)

# 绘制上面 2x2 的实线网格
for i in range(3):  # 仅绘制前两列
    ax.axvline(x[i], ymin=0.5, ymax=1, color='black', linestyle='-', linewidth=1.5)
for j in range(3):
    if j == 0:
        # j = 0 时绘制虚线
        ax.axhline(y[j], xmin=0, xmax=1, color='black', linestyle='-', linewidth=1.5)
    else:
        # j = 1 和 j = 2 时绘制实线
        ax.axhline(y[j], xmin=0, xmax=1, color='black', linestyle='-', linewidth=1.5)

# 绘制下面 2x1 的虚线网格
for i in range(3):  # 仅绘制前两列
    ax.axvline(x[i], ymin=0, ymax=0.5, color='black', linestyle='-', linewidth=1.5)
ax.axhline(y[0]/2, xmin=0, xmax=1, color='black', linestyle='-', linewidth=1.5)

# 绘制网格点
ax.plot(X, Y, 'ko', markersize=8)  # 绘制所有网格点

# 标记点 P
ax.text(1.06, 0.93, 'P', fontsize=12, color='black', ha='center', va='center')  # 中点
ax.text(1.22, 1.12, 'Q', fontsize=12, color='black', ha='center', va='center')  # 中点

# 绘制圆
circle = patches.Circle((2.3, 2.5), radius=1.8, edgecolor='black', facecolor='none', linewidth=1.5)
ax.add_patch(circle)

# 绘制从圆心到点 P 的射线
center = (2.3, 2.5)  # 圆心
point_P = (1, 1)  # 点 P
# 斜率 k
k = (point_P[1] - center[1]) / (point_P[0] - center[0])
# 截距 b
b = center[1] - k * center[0]

# 当 y = 0 时，求 x
x_intersect = -b / k
intersection_point = (x_intersect, 0)  # 射线与 y = 0 的交点

# 绘制射线
ax.plot([center[0], intersection_point[0]], [center[1], intersection_point[1]], color='black', linestyle='-', linewidth=1.5)

# 标注交点 (x, 0)
ax.plot(intersection_point[0], intersection_point[1], 'ko', markersize=8)  # 绘制交点
ax.text(intersection_point[0] + 0.1, intersection_point[1], 'T', fontsize=12, color='black', ha='left', va='bottom')

# 标注点 (1.776, 1.895)
point = (1.12, 1.14)
ax.plot(point[0], point[1], 'ko', markersize=8)  # 绘制点


# 去掉坐标轴刻度和边框
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

# 保存图片
plt.savefig('2x3_grid_with_ray.png', dpi=300)

# 显示图形
plt.show()