import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import seaborn as sns

def plot_2d_heatmap(filename, output_dir="../plot"):
    """
    从CSV文件读取数据并绘制heatmap
    
    参数:
        filename: 输入的CSV文件名
        output_dir: 输出目录
    """
    # 1. 读取CSV文件，确保使用逗号作为分隔符
    try:
        data = np.loadtxt(filename, delimiter=',')
    except ValueError as e:
        print(f"Error reading {filename}: {e}")
        return
    
    # 2. 检查数据形状，确保是二维数据
    if len(data.shape) != 2:
        print(f"Skipping irregular grid file: {filename}")
        return
    
    # 3. 绘制heatmap
    plt.figure(figsize=(10, 8))
    
    # 从文件名中提取信息
    basename = os.path.basename(filename)
    n = basename.split('_')[2][1:]  # 网格大小
    bc_type = basename.split('_')[3].split('.')[0]  # 边界条件类型
    
    # 使用 seaborn 绘制 heatmap
    sns.heatmap(data, 
                xticklabels=False,
                yticklabels=False,
                cmap='viridis',
                cbar_kws={'label': 'Solution Value'})
    
    # 4. 设置标题和坐标轴
    plt.title(f'2D Solution Heatmap (Grid Size: {n}x{n}, BC={bc_type})')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    
    # 5. 保存图形
    os.makedirs(output_dir, exist_ok=True)
    output_png = os.path.join(output_dir, basename.replace('.csv', '.png'))
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Heatmap saved to {output_png}")

# 处理所有二维规则网格文件
def process_all_2d_files(input_dir="../data"):
    files = glob.glob(os.path.join(input_dir, "solution_dim2_*.csv"))
    for file in files:
        if "irregular" not in file:  # 忽略含有 "irregular" 的文件
            plot_2d_heatmap(file)

if __name__ == "__main__":
    process_all_2d_files()
