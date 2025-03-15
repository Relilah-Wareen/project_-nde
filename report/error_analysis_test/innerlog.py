import numpy as np
import matplotlib.pyplot as plt

# 网格大小和对应的 L2 误差
grid_sizes = [8, 16, 32, 64]
l2_errors = [0.0002846379279, 7.091753574e-05, 1.824047268e-05, 6.075694405e-06]

# 转换为对数尺度
log_grid_sizes = np.log10(grid_sizes)
log_l2_errors = np.log10(l2_errors)

# 进行线性回归
slope, intercept = np.polyfit(log_grid_sizes, log_l2_errors, 1)

# 绘制对数图
plt.figure(figsize=(8, 6))
plt.plot(log_grid_sizes, log_l2_errors, 'bo-', label='Real L2 Error')

# 绘制线性回归线
regression_line = slope * np.array(log_grid_sizes) + intercept
plt.plot(log_grid_sizes, regression_line, 'r--', label=f'Linear Fit (slope = {slope:.2f})')

plt.xlabel('log10(Grid Size)')
plt.ylabel('log10(L2 Error)')
plt.title('Log-Log Plot of Inner Point with Linear Regression')
plt.grid(True, which="both", ls="--")

# 显示图例
plt.legend()

# 显示图像
plt.show()