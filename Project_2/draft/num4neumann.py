import numpy as np

# 读取C++生成的矩阵和向量
A = np.loadtxt('matrix_A.csv', delimiter=',')
F = np.loadtxt('vector_F.csv')

# --- 处理Neumann边界导致的奇异问题 ---
# 方法1：伪逆求解（适用于最小范数解）
X = np.linalg.pinv(A) @ F

# 方法2：最小二乘法（推荐）
# X = np.linalg.lstsq(A, F, rcond=None)[0]

# 方法3：添加正则化（固定一个点，消除奇异性）
# A[0, :] = 0  # 将第一行设为Dirichlet条件（如A[0,0]=1, A[0,1:]=0）
# A[0, 0] = 1
# F[0] = 0     # 固定解X[0]=0
# X = np.linalg.solve(A, F)

# 保存结果
np.savetxt('solution_X.csv', X, fmt='%.8f')

# 打印结果（注意解可能相差一个常数）
print("Solution X (Neumann边界下解可能相差一个常数):\n", X)