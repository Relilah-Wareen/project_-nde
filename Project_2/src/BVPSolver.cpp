#include "BVPSolver.hpp"
#include <stdexcept>
#include <cmath>
double IrF(double x) {
    return 1.0/16*sin(M_PI*x);
}

template <int Dim>
Vector BVPSolver<Dim>::plus_vec(const Vector& v1, const Vector& v2) {
    // 检查维度是否匹配
    if (v1.size() != v2.size()) {
        throw std::runtime_error("Vector dimensions mismatch in plus_vec");
    }
    
    Vector result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

template <int Dim>
BVPSolver<Dim>::BVPSolver(const Operator<Dim>& restrict, 
        const Operator<Dim>& prolong,
        const int& n,
        const BoundCon& bc,  // 直接传入BoundCon对象
        const Function& f,
        const int& nu1,
        const int& nu2,
        const double& epsilon)
    : restriction(restrict),
      prolongation(prolong),
      boundary_cond(bc),
      func(f),
      GridNum(n),
      nu1(nu1),
      nu2(nu2),
      err(epsilon)
{
  cnt=0;
  if(Dim != bc.getDim()) throw std::runtime_error("Solver and BoundCon dimensions mismatch.");
}

template <int Dim>
void BVPSolver<Dim>::BuildLinearSystem() {
  if constexpr (Dim == 1) {
      for(int n=GridNum; n>=2;n=n/2){
        Matrix curr_A(n+1,n+1);
        double h = 1.0 / n;
        double h2= h*h;
        for(int i=1;i<n;i++){
          curr_A.set(i,i-1, -1.0/h2);
          curr_A.set(i,i, 2.0/h2);
          curr_A.set(i,i+1, -1.0/h2);
        }
        if(boundary_cond.getBCType()==BCType::Dirichlet){
            curr_A.set(0,0,1);
            curr_A.set(n,n,1);
        }
        else if(boundary_cond.getBCType()==BCType::Neumann){
            curr_A.set(0,0, -3.0/2/h);
            curr_A.set(0,1, 2.0/h);
            curr_A.set(0,2, -1.0/2/h);
            curr_A.set(n,n-2,1.0/2/h);
            curr_A.set(n,n-1,-2.0/h);
            curr_A.set(n,n,3.0/2/h);
        }
        else if(boundary_cond.getMixed1DType()==Mixed1DType::DN){
            curr_A.set(0,0,1);
            curr_A.set(n,n-2,1.0/2/h);
            curr_A.set(n,n-1,-2.0/h);
            curr_A.set(n,n,3.0/2/h);
        }
        else if(boundary_cond.getMixed1DType()==Mixed1DType::ND){
            curr_A.set(0,0, -3.0/2/h);
            curr_A.set(0,1, 2.0/h);
            curr_A.set(0,2, -1.0/2/h);
            curr_A.set(n,n,1);
        }
        A[n]=curr_A;
      }
      int n = GridNum;
      double h = 1.0 / n;
      Vector curr_F;
      curr_F.resize(n+1,0.0);
      for(int i=0;i<curr_F.size()-1;i++){
        curr_F[i] = func(i*h);
      }
      if(boundary_cond.getBCType()==BCType::Dirichlet){
        curr_F[0]=boundary_cond.getLeft1D();
        curr_F[n]=boundary_cond.getRight1D();
    }
    else if(boundary_cond.getBCType()==BCType::Neumann){
        curr_F[0] = boundary_cond.getLeft1D();
        curr_F[n] = boundary_cond.getRight1D();
    }
    else if(boundary_cond.getMixed1DType()==Mixed1DType::DN){
        curr_F[0]=boundary_cond.getLeft1D();
        curr_F[n] = boundary_cond.getRight1D();
    }
    else if(boundary_cond.getMixed1DType()==Mixed1DType::ND){
        curr_F[0] = boundary_cond.getLeft1D();
        curr_F[n]=boundary_cond.getRight1D();
    }
    F=curr_F;

      
  } 
  else if constexpr (Dim == 2) {
    // 2D case
    if(boundary_cond.IsIrregular()) {
        BuildForIrr();
        return;
    }
    for (int n = GridNum; n >= 2; n = n / 2) {
        // 定义二维矩阵 A
        Matrix curr_A((n+1)*(n+1), (n+1)*(n+1));
        double h = 1.0 / n;
        double h2 = h * h;
        
        // 使用五点差分格式填充矩阵 A
        for (int j = 1; j < n; ++j) {
            for (int i = 1; i < n; ++i) {
                // 计算二维网格的索引 idx_2d(n, i, j)
                int idx = idx_2d(n, i, j);
                // 对角线项
                curr_A.set(idx, idx, 4.0 / h2);
                
                // 上下左右邻居项
                curr_A.set(idx, idx_2d(n, i-1, j), -1.0 / h2);  // 左
                curr_A.set(idx, idx_2d(n, i+1, j), -1.0 / h2);  // 右
                curr_A.set(idx, idx_2d(n, i, j-1), -1.0 / h2);  // 下
                curr_A.set(idx, idx_2d(n, i, j+1), -1.0 / h2);  // 上
            }
        }

        // 边界条件处理
        for (int i = 1; i < n; ++i) {
            // 处理左边界
            if (boundary_cond.getLeft2DType() == BCType::Dirichlet) {
                curr_A.set(idx_2d(n, 0, i), idx_2d(n, 0, i), 1.0);
            } else {
                // 法向导数近似
                curr_A.set(idx_2d(n, 0, i), idx_2d(n, 0, i), 1.5 / h);
                curr_A.set(idx_2d(n, 0, i), idx_2d(n, 1, i), -2.0 / h);
                curr_A.set(idx_2d(n, 0, i), idx_2d(n, 2, i), 0.5 / h);
            }
            
            // 处理右边界
            if (boundary_cond.getRight2DType() == BCType::Dirichlet) {
                curr_A.set(idx_2d(n, n, i), idx_2d(n, n, i), 1.0);
            } else {
                // 法向导数近似
                curr_A.set(idx_2d(n, n, i), idx_2d(n, n, i), 1.5 / h);
                curr_A.set(idx_2d(n, n, i), idx_2d(n, n-1, i), -2.0 / h);
                curr_A.set(idx_2d(n, n, i), idx_2d(n, n-2, i), 0.5 / h);
            }

            // 处理下边界
            if (boundary_cond.getBottom2DType() == BCType::Dirichlet) { 
                curr_A.set(idx_2d(n, i, 0), idx_2d(n, i, 0), 1.0);
            } else {
                // 法向导数近似
                curr_A.set(idx_2d(n, i, 0), idx_2d(n, i, 0), 1.5 / h);
                curr_A.set(idx_2d(n, i, 0), idx_2d(n, i, 1), -2.0 / h);
                curr_A.set(idx_2d(n, i, 0), idx_2d(n, i, 2), 0.5 / h);
            }

            // 处理上边界
            if (boundary_cond.getTop2DType() == BCType::Dirichlet) {
                curr_A.set(idx_2d(n, i, n), idx_2d(n, i, n), 1.0);
            } else {
                // 法向导数近似
                curr_A.set(idx_2d(n, i, n), idx_2d(n, i, n), 1.5 / h);
                curr_A.set(idx_2d(n, i, n), idx_2d(n, i, n-1), -2.0 / h);
                curr_A.set(idx_2d(n, i, n), idx_2d(n, i, n-2), 0.5 / h);
            }
        }
        if(boundary_cond.getBottom2DType() == BCType::Dirichlet){
            curr_A.set(idx_2d(n, 0, 0), idx_2d(n, 0, 0), 1.0);
            curr_A.set(idx_2d(n, n, 0), idx_2d(n, n, 0), 1.0);
        }
        else {
            curr_A.set(idx_2d(n, 0, 0), idx_2d(n, 0, 0), 1.5 / h);
            curr_A.set(idx_2d(n, 0, 0), idx_2d(n, 0, 1), -2.0 / h);
            curr_A.set(idx_2d(n, 0, 0), idx_2d(n, 0, 2), 0.5 / h);
            curr_A.set(idx_2d(n, n, 0), idx_2d(n, n, 0), 1.5 / h);
            curr_A.set(idx_2d(n, n, 0), idx_2d(n, n, 1), -2.0 / h);
            curr_A.set(idx_2d(n, n, 0), idx_2d(n, n, 2), 0.5 / h);
        }
        if(boundary_cond.getTop2DType() == BCType::Dirichlet){
            curr_A.set(idx_2d(n, 0, n), idx_2d(n, 0, n), 1.0);
            curr_A.set(idx_2d(n, n, n), idx_2d(n, n, n), 1.0);
        }
        else {
            curr_A.set(idx_2d(n, 0, n), idx_2d(n, 0, n), 1.5 / h);
            curr_A.set(idx_2d(n, 0, n), idx_2d(n, 0, n-1), -2.0 / h);
            curr_A.set(idx_2d(n, 0, n), idx_2d(n, 0, n-2), 0.5 / h);
            curr_A.set(idx_2d(n, n, n), idx_2d(n, n, n), 1.5 / h);
            curr_A.set(idx_2d(n, n, n), idx_2d(n, n, n-1), -2.0 / h);
            curr_A.set(idx_2d(n, n, n), idx_2d(n, n, n-2), 0.5 / h);
        }
        A[n] = curr_A;

    }
    
    // 右端项 F 的构建
    int n = GridNum;
    Vector curr_F((n+1)*(n+1), 0.0);
    double h = 1.0/n;
    for (int j = 1; j < n; ++j) {
        for (int i = 1; i < n; ++i) {
            // 对内部点赋值
            curr_F[idx_2d(n, i, j)] = func(i * h, j * h);
        }
    }
    
    // 处理边界点的右端项
    for (int i = 0; i < n+1; ++i) {
    curr_F[idx_2d(n, i, 0)] = boundary_cond.getBound2D(i * h,0);
    curr_F[idx_2d(n, i, n)] = boundary_cond.getBound2D(i * h,1);
    curr_F[idx_2d(n, 0, i)] = boundary_cond.getBound2D(0,i * h);
    curr_F[idx_2d(n, n, i)] = boundary_cond.getBound2D(1,i * h);   
    }
    F = curr_F; 
}
  else {
    throw std::runtime_error("Unexpected Dim.");
   }
}

template <int Dim>
void BVPSolver<Dim>::BuildForIrr() {
    for (int n = GridNum; n >= 2; n = n / 2) {
        // 定义二维矩阵 A
        Matrix curr_A((n+1)*(n+1), (n+1)*(n+1));
        double h = 1.0 / n;
        double h2 = h * h;
        
        // 对内部点，逐点计算其对应参数
        for (int j = 1; j < n; ++j) {
            for (int i = 1; i < n; ++i) {
                // 扭曲网格 (i*h, IrF(i*h) + j*h*(1-IrF(i*h)) );
                int idx = idx_2d(n, i, j);
                std::vector<std::vector<double>> tylor(6, std::vector<double>(6,0.0));
                Vector tgt = {0,0,0,-1,-1,0};
                Vector t(6,0.0);
                for(int k=0;k<5;k++) { // 这里的顺序也是先中心，后左右下上
                    if(k==0){
                        tylor[0][0]=1;
                    }
                    if(k==1){ // i-1, j
                        tylor[0][k]=1;
                        tylor[1][k]=-h;
                        tylor[2][k]=(1-j*h)*(IrF((i-1)*h) - IrF(i*h));
                        tylor[3][k]=0.5*tylor[1][k]*tylor[1][k];
                        tylor[4][k]=0.5*tylor[2][k]*tylor[2][k];
                        tylor[5][k]=tylor[1][k]*tylor[2][k];
                    }
                    if(k==2){ // i+1,j
                        tylor[0][k]=1;
                        tylor[1][k]=h;
                        tylor[2][k]=(1-j*h)*(IrF((i+1)*h) - IrF(i*h));
                        tylor[3][k]=0.5*tylor[1][k]*tylor[1][k];
                        tylor[4][k]=0.5*tylor[2][k]*tylor[2][k];
                        tylor[5][k]=tylor[1][k]*tylor[2][k];
                    }
                    if(k==3){ // i,j-1
                        tylor[0][k]=1;
                        tylor[1][k]=0;
                        tylor[2][k]=h*(IrF(i*h)-1);
                        tylor[3][k]=0.5*tylor[1][k]*tylor[1][k];
                        tylor[4][k]=0.5*tylor[2][k]*tylor[2][k];
                        tylor[5][k]=tylor[1][k]*tylor[2][k];
                    }
                    if(k==4){ // i,j+1
                        tylor[0][k]=1;
                        tylor[1][k]=0;
                        tylor[2][k]=-h*(IrF(i*h)-1);
                        tylor[3][k]=0.5*tylor[1][k]*tylor[1][k];
                        tylor[4][k]=0.5*tylor[2][k]*tylor[2][k];
                        tylor[5][k]=tylor[1][k]*tylor[2][k];
                    }
                }
                // i-1, j-1
                tylor[0][5]=1;
                tylor[1][5]=-h;
                tylor[2][5]=(1-(j-1)*h)*IrF((i-1)*h)-h-(1-j*h)*IrF(i*h);
                tylor[3][5]=0.5*tylor[1][5]*tylor[1][5];
                tylor[4][5]=0.5*tylor[2][5]*tylor[2][5];
                tylor[5][5]=tylor[1][5]*tylor[2][5];
                t = LUSolver(tylor,tgt);
                curr_A.add(idx, idx, 0.25*t[0]);                
                curr_A.add(idx, idx_2d(n, i-1, j), 0.25*t[1]);  // 左
                curr_A.add(idx, idx_2d(n, i+1, j), 0.25*t[2]);  // 右
                curr_A.add(idx, idx_2d(n, i, j-1), 0.25*t[3]);  // 下
                curr_A.add(idx, idx_2d(n, i, j+1), 0.25*t[4]);  // 上
                curr_A.add(idx, idx_2d(n,i-1,j-1), 0.25*t[5]);
                // i+1, j-1
                tylor[0][5]=1;
                tylor[1][5]=h;
                tylor[2][5]=(1-(j-1)*h)*IrF((i+1)*h)-h-(1-j*h)*IrF(i*h);
                tylor[3][5]=0.5*tylor[1][5]*tylor[1][5];
                tylor[4][5]=0.5*tylor[2][5]*tylor[2][5];
                tylor[5][5]=tylor[1][5]*tylor[2][5];
                t = LUSolver(tylor,tgt);
                curr_A.add(idx, idx, 0.25*t[0]);                
                curr_A.add(idx, idx_2d(n, i-1, j), 0.25*t[1]);  // 左
                curr_A.add(idx, idx_2d(n, i+1, j), 0.25*t[2]);  // 右
                curr_A.add(idx, idx_2d(n, i, j-1), 0.25*t[3]);  // 下
                curr_A.add(idx, idx_2d(n, i, j+1), 0.25*t[4]);  // 上
                curr_A.add(idx, idx_2d(n,i+1,j-1), 0.25*t[5]);
                // i-1, j+1
                tylor[0][5]=1;
                tylor[1][5]=-h;
                tylor[2][5]=(1-(j+1)*h)*IrF((i-1)*h)+h-(1-j*h)*IrF(i*h);
                tylor[3][5]=0.5*tylor[1][5]*tylor[1][5];
                tylor[4][5]=0.5*tylor[2][5]*tylor[2][5];
                tylor[5][5]=tylor[1][5]*tylor[2][5];
                t = LUSolver(tylor,tgt);
                curr_A.add(idx, idx, 0.25*t[0]);                
                curr_A.add(idx, idx_2d(n, i-1, j), 0.25*t[1]);  // 左
                curr_A.add(idx, idx_2d(n, i+1, j), 0.25*t[2]);  // 右
                curr_A.add(idx, idx_2d(n, i, j-1), 0.25*t[3]);  // 下
                curr_A.add(idx, idx_2d(n, i, j+1), 0.25*t[4]);  // 上
                curr_A.add(idx, idx_2d(n,i-1,j+1), 0.25*t[5]);
                // i+1, j+1
                tylor[0][5]=1;
                tylor[1][5]=h;
                tylor[2][5]=(1-(j+1)*h)*IrF((i+1)*h)+h-(1-j*h)*IrF(i*h);
                tylor[3][5]=0.5*tylor[1][5]*tylor[1][5];
                tylor[4][5]=0.5*tylor[2][5]*tylor[2][5];
                tylor[5][5]=tylor[1][5]*tylor[2][5];
                t = LUSolver(tylor,tgt);
                curr_A.add(idx, idx, 0.25*t[0]);                
                curr_A.add(idx, idx_2d(n, i-1, j), 0.25*t[1]);  // 左
                curr_A.add(idx, idx_2d(n, i+1, j), 0.25*t[2]);  // 右
                curr_A.add(idx, idx_2d(n, i, j-1), 0.25*t[3]);  // 下
                curr_A.add(idx, idx_2d(n, i, j+1), 0.25*t[4]);  // 上
                curr_A.add(idx, idx_2d(n,i+1,j+1), 0.25*t[5]);
            }
        }

        // 边界条件处理
        for (int i = 1; i < n; ++i) {
            // 处理左边界
            if (boundary_cond.getLeft2DType() == BCType::Dirichlet) {
                curr_A.set(idx_2d(n, 0, i), idx_2d(n, 0, i), 1.0);
            } else {
                throw std::runtime_error("Irregular only support Dirichlet bound.");
            }
            
            // 处理右边界
            if (boundary_cond.getRight2DType() == BCType::Dirichlet) {
                curr_A.set(idx_2d(n, n, i), idx_2d(n, n, i), 1.0);
            } else {
                throw std::runtime_error("Irregular only support Dirichlet bound.");
            }

            // 处理下边界
            if (boundary_cond.getBottom2DType() == BCType::Dirichlet) { 
                curr_A.set(idx_2d(n, i, 0), idx_2d(n, i, 0), 1.0);
            } else {
                throw std::runtime_error("Irregular only support Dirichlet bound.");
            }

            // 处理上边界
            if (boundary_cond.getTop2DType() == BCType::Dirichlet) {
                curr_A.set(idx_2d(n, i, n), idx_2d(n, i, n), 1.0);
            } else {
                throw std::runtime_error("Irregular only support Dirichlet bound.");
            }
        }
        if(boundary_cond.getLeft2DType() == BCType::Dirichlet){
            curr_A.set(idx_2d(n, 0, 0), idx_2d(n, 0, 0), 1.0);
            curr_A.set(idx_2d(n, 0, n), idx_2d(n, 0, n), 1.0);
        }
        else {
            throw std::runtime_error("Irregular only support Dirichlet bound.");
        }
        if(boundary_cond.getRight2DType() == BCType::Dirichlet){
            curr_A.set(idx_2d(n, n, 0), idx_2d(n, n, 0), 1.0);
            curr_A.set(idx_2d(n, n, n), idx_2d(n, n, n), 1.0);
        }
        else {
            throw std::runtime_error("Irregular only support Dirichlet bound.");
        }
        A[n] = curr_A;
    }
    
    // 右端项 F 的构建
    int n = GridNum;
    Vector curr_F((n+1)*(n+1), 0.0);
    double h = 1.0/n;
    for (int j = 1; j < n; ++j) {
        for (int i = 1; i < n; ++i) {
            // 对内部点赋值
            curr_F[idx_2d(n, i, j)] = func(i*h, IrF(i*h) + j*h*(1-IrF(i*h)) );
        }
    }
    
    // 处理边界点的右端项
    for (int i = 0; i < n+1; ++i) {
    curr_F[idx_2d(n, i, 0)] = boundary_cond.getBound2D(i * h,IrF(i*h));
    curr_F[idx_2d(n, i, n)] = boundary_cond.getBound2D(i * h,1);
    curr_F[idx_2d(n, 0, i)] = boundary_cond.getBound2D(0,i * h);
    curr_F[idx_2d(n, n, i)] = boundary_cond.getBound2D(1,i * h);   
    }
    F = curr_F;
}

template <int Dim>
Vector BVPSolver<Dim>::LUSolver(const std::vector<std::vector<double>>& A, const Vector& b) const {
    int n = A.size();
    if (n == 0 || A[0].size() != n || b.size() != n) {
        throw std::runtime_error("Invalid matrix or vector dimensions in LUSolver");
    }

    // 创建增广矩阵的拷贝
    std::vector<std::vector<double>> LU = A;
    Vector x(n, 0.0);

    // LU分解
    for (int k = 0; k < n; ++k) {
        // 计算U的行
        for (int j = k; j < n; ++j) {
            double sum = 0.0;
            for (int p = 0; p < k; ++p) {
                sum += LU[k][p] * LU[p][j];
            }
            LU[k][j] = A[k][j] - sum;
        }

        // 计算L的列
        for (int i = k + 1; i < n; ++i) {
            double sum = 0.0;
            for (int p = 0; p < k; ++p) {
                sum += LU[i][p] * LU[p][k];
            }
            LU[i][k] = (A[i][k] - sum) / LU[k][k];
        }
    }

    // 前向替换（Ly = b）
    Vector y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int p = 0; p < i; ++p) {
            sum += LU[i][p] * y[p];
        }
        y[i] = b[i] - sum;
    }

    // 后向替换（Ux = y）
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int p = i + 1; p < n; ++p) {
            sum += LU[i][p] * x[p];
        }
        x[i] = (y[i] - sum) / LU[i][i];
    }

    return x;
}

template <int Dim>
void BVPSolver<Dim>::SolveByLU() {
    auto it_A = A.find(GridNum);
    if (it_A == A.end()) {
        throw std::runtime_error("Matrix or vector not built for current grid");
    }

    // 将稀疏矩阵转换为稠密矩阵
    const Matrix& sparseA = it_A->second;
    int n = sparseA.getNumRows();
    std::vector<std::vector<double>> denseA(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            denseA[i][j] = sparseA.get(i, j);
        }
    }

    // 调用LU求解器
    result = LUSolver(denseA, F);
}

template <int Dim>
void BVPSolver<Dim>::SolveByRelaxOnly(double omega) {
    auto it_A = A.find(GridNum);
    if (it_A == A.end()) {
        throw std::runtime_error("Matrix or vector not built for current grid");
    }
    
    const Matrix& mat = it_A->second;
    const Vector& curr_F = F;
    Vector x(mat.getNumCols(), 0.0); // 初始猜测设为0
    
    const int max_iter = 100000;
    bool converged = false;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        x = relax(mat, x, curr_F, omega);
        
        // 计算残差
        Vector r = residual(mat, x, curr_F);
        
        // 计算残差的无穷范数
        double res_norm = 0.0;
        for (double val : r) {
            res_norm = std::max(res_norm, std::abs(val));
        }
        
        // 检查收敛
        if (res_norm < err) {
            converged = true;
            // std::cout << "Converged after " << iter+1 << " iterations" << std::endl;
            break;
        }
    }
    
    if (!converged) {
        std::cout << "Warning: Reached maximum iterations without convergence" << std::endl;
    }
    
    result = x;
}


template <int Dim>
void BVPSolver<Dim>::SolveByVC() {
    if constexpr (Dim == 1  || Dim == 2) {
        auto it_A = A.find(GridNum);
        if (it_A == A.end()) {
            throw std::runtime_error("Matrix or vector not built for current grid");
        }
        
        const Matrix& mat = it_A->second;
        const Vector& curr_F = F;
        Vector x = getZeroVector(GridNum);
        
        const int max_iter = 1000000;
        bool converged = false;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            x = VC(GridNum, x, curr_F);
            
            // 计算残差
            Vector r = residual(mat, x, curr_F);
            
            // 计算残差的无穷范数
            double res_norm = 0.0;
            for (double val : r) {
                res_norm = std::max(res_norm, std::abs(val));
            }
            
            // 检查收敛
            if (res_norm < err) {
                converged = true;
                // std::cout << "V-cycle converged after " << iter+1 << " iterations" << std::endl;
                break;
            }
        }
        
        if (!converged) {
            std::cout << "Warning: V-cycle reached maximum iterations without convergence" << std::endl;
        }
        
        result = x;
    } else {
        throw std::runtime_error("Unexpected Dim for Vcycle.");
    }
}

template <int Dim>
void BVPSolver<Dim>::SolveByFMG() {
    auto it_A = A.find(GridNum);
    if (it_A == A.end()) {
        throw std::runtime_error("Matrix or vector not built for current grid");
    }
    
    const Matrix& mat = it_A->second;
    const Vector& curr_F = F;
    Vector x = getZeroVector(GridNum);
    Vector f = curr_F;
    const int max_iter = 1000000;
    bool converged = false;

    for (int iter = 0; iter < max_iter; ++iter) {
        x = plus_vec(x , FMG(GridNum, f));
        
        // 计算残差
        Vector r = residual(mat, x, curr_F);
        
        // 计算残差的无穷范数
        double res_norm = 0.0;
        for (double val : r) {
            res_norm = std::max(res_norm, std::abs(val));
        }
        
        // 检查收敛
        if (res_norm < err) {
            converged = true;
            // std::cout << "FMG converged after " << iter+1 << " iterations" << std::endl;
            break;
        }
        
        // 如果不收敛，更新右端项为残差
        f = r;
    }
    
    if (!converged) {
        std::cout << "Warning: FMG reached maximum iterations without convergence" << std::endl;
    }

    result = x;

}


template <int Dim>
const Vector& BVPSolver<Dim>::getResult() const {
    return result;
}

template <int Dim>
void BVPSolver<Dim>::printA() const {
    if constexpr (Dim == 1) {
        auto it = A.find(GridNum);
        if (it != A.end()) {
            std::cout << "Matrix A (GridNum = " << GridNum << "):" << std::endl;
            it->second.print();
        } else {
            std::cout << "Matrix A for GridNum " << GridNum << " not found." << std::endl;
        }
    } else {
        throw std::runtime_error("Only Dim=1 is implemented");
    }
}

template <int Dim>
void BVPSolver<Dim>::printF() const {
    for (size_t i = 0; i < F.size(); ++i) {
        std::cout << F[i] << std::endl;
    }
 
}

template <int Dim>
Vector BVPSolver<Dim>::relax(const Matrix& A, const Vector& v, const Vector& f, double omega) {
    if (A.getNumRows() != A.getNumCols() || 
        A.getNumRows() != static_cast<int>(v.size()) || 
        A.getNumRows() != static_cast<int>(f.size())) {
        throw std::invalid_argument("Matrix and vector dimensions mismatch");
    }

    Vector x_new = v; // 复制初始猜测
    const int n = A.getNumRows();

    for (int i = 0; i < n; ++i) {
        double sigma = 0.0;
        double diag = 0.0;
        
        // 遍历矩阵第i行的非零元素
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k) {
            const int j = A.columns[k];
            if (j != i) {
                sigma += A.values[k] * x_new[j]; // 使用已更新的x_new
            } else {
                diag = A.values[k];
            }
        }
        
        if (diag == 0.0) {
            throw std::runtime_error("Zero diagonal element encountered");
        }
        
        // Gauss-Seidel with relaxation
        x_new[i] = (1.0 - omega) * x_new[i] + omega * (f[i] - sigma) / diag;
    }
    cnt++;
    return x_new;
}

template <int Dim>
Vector BVPSolver<Dim>::residual(const Matrix& A, const Vector& v, const Vector& f) const {
    if (A.getNumRows() != static_cast<int>(f.size()) || 
        A.getNumCols() != static_cast<int>(v.size())) {
        throw std::invalid_argument("Matrix and vector dimensions mismatch in residual calculation");
    }
    
    Vector r = f;
    Vector Av = A.multiply(v);
    
    for (size_t i = 0; i < r.size(); ++i) {
        r[i] -= Av[i];
    }
    
    return r;
}


template <int Dim>
void BVPSolver<Dim>::printResult() const {
    if constexpr (Dim == 1) {
        std::cout << "Solution result (GridNum = " << GridNum << "):" << std::endl;
        for (size_t i = 0; i < result.size(); ++i) {
            std::cout << i << "\t" << result[i] << std::endl;
        }
    } else if constexpr (Dim == 2) {
        std::cout << "Solution result (GridNum = " << GridNum << "):" << std::endl;
        for (int j = GridNum; j >= 0; --j) {
            for (int i = 0; i <= GridNum; ++i) {
                std::cout << result[idx_2d(GridNum, i, j)] << "\t";
            }
            std::cout << std::endl;
        }
    }
    else {
        
    }
}

template <int Dim>
void BVPSolver<Dim>::printResultToCSV(const std::string& filename) const {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::runtime_error("Failed to open file for writing");
    }

    if constexpr (Dim == 1) {
        // 一维输出：每行一个点 (x, u(x))
        for (size_t i = 0; i < result.size(); ++i) {
            double x = static_cast<double>(i) / GridNum;
            outFile << x << "," << result[i] << std::endl;
        }
    } else if constexpr (Dim == 2) {
        // 二维输出：按几何网格排列
        if(boundary_cond.IsIrregular()){
            double h = 1.0 / GridNum;
            for (int j = GridNum; j >= 0; --j) {
                for (int i = 0; i <= GridNum; ++i) {
                    outFile << i*h << ","<< IrF(i*h) + j*h*(1-IrF(i*h)) << "," << result[idx_2d(GridNum, i, j)] << std::endl;
                }
            }
        }
        else {
            for (int j = GridNum; j >= 0; --j) {
            for (int i = 0; i <= GridNum; ++i) {
                outFile << result[idx_2d(GridNum, i, j)];
                if (i != GridNum) outFile << ",";
            }
            outFile << std::endl;
        }}
    }
    outFile.close();
}


template<int Dim>
int BVPSolver<Dim>::getCount() const {
  return cnt;
}

template <int Dim>
Vector BVPSolver<Dim>::VC(const int& n, Vector v, const Vector& f) {
    double omega=1;
    if(Dim == 1) omega = 2.0/3;
    else if(Dim == 2) omega = 4.0/5;
    else {
        throw std::runtime_error("Unexpected Dim.");
    }
    // 步骤1: 松弛nu1次
    auto it = A.find(n);
    if (it == A.end()) {
        throw std::runtime_error("Matrix not built for current grid");
    }
    const Matrix& mat = it->second;
    
    for (int i = 0; i < nu1; ++i) {
        v = relax(mat, v, f, omega);
    }
    
    // 步骤2: 检查是否是最粗网格
    if (n <= 2) {
        // 步骤4: 直接进行nu2次松弛
        for (int i = 0; i < nu2; ++i) {
            v = relax(mat, v, f, omega);
        }
        return v;
    }
    
    // 计算残差并限制到粗网格
    Vector res = residual(mat, v, f);
    Vector res_coarse = restriction(res); // 使用限制算子

    Vector v_coarse = getZeroVector(n/2); // 创建粗网格零向量
    v_coarse = VC(n/2, v_coarse, res_coarse);
    
    // 步骤3: 延长并修正
    Vector correction = prolongation(v_coarse);
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] += correction[i];
    }
    
    // 步骤4: 后松弛nu2次
    for (int i = 0; i < nu2; ++i) {
        v = relax(mat, v, f, omega);
    }
    
    return v;

}

template <int Dim>
Vector BVPSolver<Dim>::FMG(const int& n, const Vector& f) {
    Vector v;
    
    // 检查是否是最粗网格
    if (n <= 2) {
        v = getZeroVector(n);
    } else {
    v = prolongation(FMG(n/2, restriction(f)));
    }
    
    // 执行V-cycle
    return VC(n, v, f);

}

template <int Dim>
Vector BVPSolver<Dim>::getZeroVector(int n) const {
  if constexpr (Dim == 1) {
      return Vector(n+1, 0.0);
  } 
  else if constexpr (Dim == 2){
    return Vector((n+1)*(n+1), 0.0);
  }
  else {
      throw std::runtime_error("Only Dim=1 is implemented");
  }
}

template <int Dim>
void BVPSolver<Dim>::saveMatrixToCSV(const std::string& filename) const {
    auto it = A.find(GridNum);
    if (it == A.end()) {
        throw std::runtime_error("Matrix not built for current grid");
    }
    
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::runtime_error("Failed to open file for writing");
    }
    
    const Matrix& mat = it->second;
    const int n = mat.getNumRows();
    
    // 写入CSV头(可选)
    // outFile << "i,j,value" << std::endl;
    
    // 写入矩阵数据
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // 使用逗号分隔，每行一个元素
            outFile << mat.get(i, j);
            if (j != n - 1) {
                outFile << ",";
            }
        }
        outFile << std::endl;
    }
    
    outFile.close();
}

template <int Dim>
void BVPSolver<Dim>::saveVectorToCSV(const std::string& filename) const {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw std::runtime_error("Failed to open file for writing");
    }
    
    const Vector& vec = F;
    
    // 写入向量数据，每行一个元素
    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i] << std::endl;
    }
    
    outFile.close();
}

template <int Dim>
void BVPSolver<Dim>::normalizeNeumannSolution() {
    if (result.empty()) {
        throw std::runtime_error("Result vector is empty. Solve the problem first.");
    }

    // 计算平移量：目标是将 result[0] 设为 1，因此平移量为 (1 - result[0])
    double shift = 1.0 - result[0];

    // 应用平移
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] += shift;
    }
}

template class BVPSolver<1>;
template class BVPSolver<2>;  // 暂不实例化Dim=2