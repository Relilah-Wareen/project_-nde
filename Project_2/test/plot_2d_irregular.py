import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import numpy as np

def u_true_2d(x, y):
    """2D 真实解函数: exp(sin(x)+y)"""
    return np.exp(np.sin(x) + y)

def plot_2d_irregular(filename, output_dir="../plot"):
    """
    绘制二维不规则网格的散点图
    
    参数:
        filename: 输入的CSV文件名
        output_dir: 输出目录
    """
    try:
        # 1. 读取数据，指定分隔符为逗号
        df = pd.read_csv(filename, header=None, names=['x', 'y', 'u'], delimiter=',')
    except pd.errors.ParserError as e:
        print(f"Error reading {filename}: {e}")
        return
    
    # 2. 检查数据格式
    if df.shape[1] != 3:
        print(f"Skipping file with invalid format: {filename}")
        return
    
    # 3. 创建图形
    plt.figure(figsize=(12, 10))
    
    # 4. 绘制散点热图
    sc = plt.scatter(df['x'], df['y'], c=df['u'], s=50, cmap='viridis')
    
    # 5. 添加颜色条
    cbar = plt.colorbar(sc)
    cbar.set_label('u(x,y)')
    
    # 6. 从文件名中提取信息
    basename = os.path.basename(filename)
    n = basename.split('_')[2][1:]  # 网格大小
    bc_type = basename.split('_')[3].split('.')[0]  # 边界条件类型
    
    # 7. 设置标题和坐标轴
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Numerical Solution on Irregular Grid (n={n}, BC={bc_type})')
    
    # 8. 保存图形
    os.makedirs(output_dir, exist_ok=True)
    output_png = os.path.join(output_dir, basename.replace('.csv', '.png'))
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Irregular grid plot saved to {output_png}")

# 处理所有二维不规则网格文件
def process_all_irregular_files(input_dir="../data"):
    files = glob.glob(os.path.join(input_dir, "*irregular*.csv"))
    for file in files:
        plot_2d_irregular(file)

if __name__ == "__main__":
    process_all_irregular_files()
