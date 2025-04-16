import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_heatmap_from_csv(filename):
    """
    从CSV文件读取数据并绘制heatmap
    
    参数:
        filename: 输入的CSV文件名
    """
    # 1. 读取CSV文件
    data = np.loadtxt(filename, delimiter=',')
    
    # 2. 创建坐标网格
    n = data.shape[1] - 1  # 网格大小 (因为点数是n+1)
    x = np.linspace(0, 1, n+1)
    y = np.linspace(0, 1, n+1)
    X, Y = np.meshgrid(x, y)
    
    # 3. 绘制heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(data, 
               xticklabels=False,  # 去除x轴坐标标签
               yticklabels=False,  # 去除y轴坐标标签
               cmap='viridis',
               cbar_kws={'label': 'Solution Value'})
    
    # 4. 设置标题和坐标轴
    plt.title('Solution Heatmap (Grid Size: {}x{})'.format(n, n))
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    
    # 5. 保存和显示
    output_png = filename.replace('.csv', '.png')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"Heatmap saved to {output_png}")

# 使用示例
plot_heatmap_from_csv("solution.csv")
