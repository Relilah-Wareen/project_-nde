import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# 网格大小和对应的边界 L2 误差
grid_sizes = [8, 16, 32, 64]
boundary_l2_errors = [0.391238705316,0.143915371632,0.0811701567764,0.01303440544]

# 转换为对数尺度
log_grid_sizes = np.log10(grid_sizes)
log_boundary_l2_errors = np.log10(boundary_l2_errors)

# 线性回归
slope, intercept, r_value, p_value, std_err = linregress(log_grid_sizes, log_boundary_l2_errors)

# 绘制对数图
plt.figure(figsize=(8, 6))
plt.plot(log_grid_sizes, log_boundary_l2_errors, 'bo-', label='Boundary L2 Error')
plt.xlabel('log10(Grid Size)')
plt.ylabel('log10(Boundary L2 Error)')
plt.title('Log-Log Plot of Boundary L2 Error with Linear Regression')
plt.grid(True, which="both", ls="--")

# 添加回归线
regression_line = slope * np.array(log_grid_sizes) + intercept
plt.plot(log_grid_sizes, regression_line, 'g--', label=f'Regression Line (Slope = {slope:.2f})')

# 显示图例
plt.legend()

# 显示图像
plt.show()

# 输出回归斜率
print(f"Regression Slope: {slope:.4f}")