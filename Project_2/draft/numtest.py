import numpy as np

# 直接读取C++生成的文件
A = np.loadtxt('matrix_A.csv', delimiter=',')
F = np.loadtxt('vector_F.csv')

# 解方程并保存结果
X = np.linalg.solve(A, F)
np.savetxt('solution_X.csv', X, fmt='%.8f')

# 打印结果
print("Solution X:\n", X)