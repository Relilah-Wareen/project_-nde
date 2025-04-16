import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取数据
df = pd.read_csv('solution.csv', header=None, names=['x', 'y', 'u'])

# 创建图形
plt.figure(figsize=(10, 8))

# 绘制散点热图（每个点的大小和颜色代表u值）
sc = plt.scatter(df['x'], df['y'], c=df['u'], s=100, cmap='viridis')

# 添加颜色条
cbar = plt.colorbar(sc)
cbar.set_label('u(x,y)')

# 标记坐标轴
plt.xlabel('x')
plt.ylabel('y')
plt.title('Numerical Solution (Direct Scatter Plot)')

plt.tight_layout()
plt.show()