#include "SparseMatrix.hpp"
#include <iostream>

void test_sparse_matrix() {
    // 创建一个 3x3 的稀疏矩阵
    SparseMatrix mat(3, 3);
    
    // 设置一些非零元素
    mat.set(0, 0, 1.0);
    mat.set(1, 1, 2.0);
    mat.set(2, 2, 3.0);
    mat.set(0, 1, 0.5);
    mat.set(1, 0, 0.5);
    
    // 打印矩阵
    std::cout << "Matrix:" << std::endl;
    mat.print();
    
    // 测试矩阵向量乘法
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y = mat.multiply(x);
    
    std::cout << "\nMatrix-vector product (A * x):" << std::endl;
    for (double val : y) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    // 测试从三重格式构造
    std::vector<int> rows = {0, 1, 2, 0, 1};
    std::vector<int> cols = {0, 1, 2, 1, 0};
    std::vector<double> vals = {1.0, 2.0, 3.0, 0.5, 0.5};
    
    SparseMatrix mat2(3, 3, rows, cols, vals);
    std::cout << "\nMatrix constructed from triplets:" << std::endl;
    mat2.print();
}

int main() {
    test_sparse_matrix();
    return 0;
}