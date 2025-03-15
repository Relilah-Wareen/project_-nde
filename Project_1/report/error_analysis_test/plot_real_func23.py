import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# 定义文件名列表
file_names = [
    'true_solution_exp.csv',
    'true_solution_linear.csv'
]

# 定义输入和输出目录
input_dir = '../../data'  # CSV 文件所在的目录
output_dir = '../error_analysis_test'  # 图像保存的目录

# 确保输出目录存在
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 遍历每个文件
for file_name in file_names:
    # 构建输入文件路径
    input_file_path = os.path.join(input_dir, file_name)
    
    # 读取CSV文件
    data = pd.read_csv(input_file_path, header=None)

    # 将数据转换为NumPy数组
    heatmap_data = data.to_numpy()

    # 绘制热力图
    plt.figure(figsize=(10, 8))
    plt.imshow(heatmap_data, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Value')
    
    # 设置标题和坐标轴标签
    plt.title(f'Heatmap of True Solution ({file_name})')
    plt.xlabel('X Axis')
    plt.ylabel('Y Axis')

    # 构建输出文件路径
    output_image_name = os.path.join(output_dir, file_name.replace('.csv', '.png'))
    # 保存图像到 ../plot 目录
    plt.savefig(output_image_name, bbox_inches='tight', dpi=300)
    plt.close() 