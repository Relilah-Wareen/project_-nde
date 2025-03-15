import numpy as np
import matplotlib.pyplot as plt

# 网格大小和对应的边界 L2 误差
grid_sizes = [8, 16, 32, 64]
boundary_l2_errors = [0.02517359827, 0.007352487249, 0.002116067993, 0.0005991869491]

# 转换为对数尺度
log_grid_sizes = np.log10(grid_sizes)
log_boundary_l2_errors = np.log10(boundary_l2_errors)

# 进行线性回归
slope, intercept = np.polyfit(log_grid_sizes, log_boundary_l2_errors, 1)

# 绘制对数图
plt.figure(figsize=(8, 6))
plt.plot(log_grid_sizes, log_boundary_l2_errors, 'bo-', label='Boundary L2 Error')

# 绘制线性回归线
regression_line = slope * np.array(log_grid_sizes) + intercept
plt.plot(log_grid_sizes, regression_line, 'r--', label=f'Linear Fit (slope = {slope:.2f})')

plt.xlabel('log10(Grid Size)')
plt.ylabel('log10(Boundary L2 Error)')
plt.title('Log-Log Plot of Neumann Square Boundary L2 Error with Linear Regression')
plt.grid(True, which="both", ls="--")

# 显示图例
plt.legend()

# 显示图像
plt.show()