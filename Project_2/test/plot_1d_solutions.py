import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def u_true_1d(x):
    """1D 真实解函数: exp(sin(x)+x)"""
    return np.exp(np.sin(x)+x)

def plot_1d_solution(filename, output_dir="../plot"):
    """
    绘制一维数值解与真实解的比较图
    
    参数:
        filename: 输入的CSV文件名
        output_dir: 输出目录
    """
    # 读取数据，指定逗号作为分隔符
    try:
        data = np.loadtxt(filename, delimiter=',')
    except ValueError as e:
        print(f"Error reading {filename}: {e}")
        return

    # 数据拆分
    x = data[:, 0]
    u_num = data[:, 1]
    
    # 计算真实解
    u_true = u_true_1d(x)
    
    # 创建图形
    plt.figure(figsize=(10, 6))
    
    # 绘制数值解和真实解
    plt.plot(x, u_num, 'b-', linewidth=2, label='Numerical Solution')
    plt.plot(x, u_true, 'r--', linewidth=2, label='True Solution')
    
    # 设置标题和标签
    basename = os.path.basename(filename)
    n = basename.split('_')[2][1:]  # 从文件名中提取网格大小
    bc_type = basename.split('_')[3].split('.')[0]  # 提取边界条件类型
    
    plt.title(f'1D Solution Comparison (n={n}, BC={bc_type})')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.legend()
    plt.grid(True)
    
    # 保存图形
    os.makedirs(output_dir, exist_ok=True)
    output_png = os.path.join(output_dir, basename.replace('.csv', '.png'))
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {output_png}")

# 处理所有一维数据文件
def process_all_1d_files(input_dir="../data"):
    files = glob.glob(os.path.join(input_dir, "solution_dim1_*.csv"))
    for file in files:
        plot_1d_solution(file)

if __name__ == "__main__":
    process_all_1d_files()
