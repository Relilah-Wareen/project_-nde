import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

# 配置路径
DATA_DIR = "../data"
PLOT_DIR = "../plot"
os.makedirs(PLOT_DIR, exist_ok=True)

def plot_phase_portrait(file_path, output_dir):
    """为单个CSV文件生成相位图"""
    try:
        # 解析文件名
        filename = Path(file_path).stem
        parts = filename.split("_")
        
        # 检查文件名格式是否合法
        if len(parts) < 4:
            print(f"文件名格式错误: {filename} (至少需要4个下划线分隔的部分)")
            return

        model_name = parts[0]
        method_name = parts[1]
        p_value = parts[2][1:] if parts[2].startswith("p") else "unknown"
        n_value = parts[3][1:] if parts[3].startswith("n") else "unknown"

        # 创建方法专属目录（确保路径安全）
        method_dir = os.path.join(output_dir, method_name)
        os.makedirs(method_dir, exist_ok=True)

        # 读取数据
        df = pd.read_csv(file_path)
        if 'u1' not in df.columns or 'u2' not in df.columns:
            print(f"CSV文件缺少u1/u2列: {filename}")
            return

        u1 = df['u1']
        u2 = df['u2']

        # 创建画布
        plt.figure(figsize=(10, 8))

        # 绘制相位轨迹
        plt.plot(u1, u2, 'b-', linewidth=1.5, alpha=0.8)
        plt.scatter(u1[0], u2[0], c='lime', s=80, edgecolor='k', label='Start')
        plt.scatter(u1.iloc[-1], u2.iloc[-1], c='r', s=80, edgecolor='k', label='End')

        # 添加标注信息
        info_text = f"{model_name}\np={p_value}, n={n_value}"
        plt.annotate(info_text, xy=(0.05, 0.95), xycoords='axes fraction',
                     fontsize=12, ha='left', va='top',
                     bbox=dict(boxstyle='round', alpha=0.2, facecolor='w'))

        # 设置图表属性
        plt.title(f'{method_name} Method - Phase Portrait', fontsize=14)
        plt.xlabel('u1', fontsize=12)
        plt.ylabel('u2', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.axis('equal')
        plt.legend(loc='upper right')

        # 保存文件
        output_path = os.path.join(method_dir, f"{filename}.png")
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Generated: {output_path}")

    except Exception as e:
        print(f"处理文件 {filename} 时发生错误: {str(e)}")

if __name__ == "__main__":
    # 遍历所有CSV文件
    for filename in os.listdir(DATA_DIR):
        if filename.endswith(".csv"):
            file_path = os.path.join(DATA_DIR, filename)
            plot_phase_portrait(file_path, PLOT_DIR)