import numpy as np

def solve_and_save(grid_num, matrix_file='matrix_A.csv', vector_file='vector_F.csv', output_file='solution_X.csv'):
    """
    解方程并按几何位置格式输出结果
    
    参数:
        grid_num: 网格划分数 (如4表示5x5个点)
        matrix_file: 矩阵文件路径
        vector_file: 向量文件路径 
        output_file: 输出文件路径
    """
    # 1. 读取数据
    A = np.loadtxt(matrix_file, delimiter=',')
    F = np.loadtxt(vector_file)
    
    # 2. 解方程
    X = np.linalg.solve(A, F)
    
    # 3. 按几何位置重新排列 (二维情况)
    if len(F) == (grid_num + 1)**2:  # 判断是否为二维问题
        # 重新组织为网格结构
        solution_grid = X.reshape((grid_num + 1, grid_num + 1))
        
        # 按从下到上的顺序输出 (y坐标递增)
        with open(output_file, 'w') as f:
            for j in range(grid_num + 1):
                row = solution_grid[:, j]  # 取一列数据 (x方向)
                f.write(','.join(f"{val:.8f}" for val in row) + '\n')
        
        print(f"Solution saved to {output_file} in grid layout")
        return solution_grid
    
    # 一维情况直接保存
    np.savetxt(output_file, X, fmt='%.8f')
    print(f"Solution saved to {output_file}")
    return X

# 使用示例 (假设是4x4网格，共5x5=25个点)
grid_size = 8  # 根据实际网格大小修改
solution = solve_and_save(grid_size)

# 可视化 (可选)
import matplotlib.pyplot as plt
if isinstance(solution, np.ndarray) and solution.ndim == 2:
    plt.imshow(solution, cmap='hot', origin='lower')
    plt.colorbar(label='Solution Value')
    plt.title('Solution Heatmap')
    plt.xlabel('X index')
    plt.ylabel('Y index')
    plt.show()