#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

#include <vector>
#include <cstddef> 
typedef std::vector<double> Vector; 

class SparseMatrix {
public:
     // 定义Vector别名
    Vector values;                      // 非零元素值
    std::vector<int> columns;           // 非零元素的列索引
    std::vector<int> row_ptr;           // 行指针
    int num_rows;                       // 矩阵行数
    int num_cols;                       // 矩阵列数
    size_t nnz;                         // 非零元素个数

    void insert(int row, int col, double value); // 辅助插入函数
    // 构造函数
    SparseMatrix() : num_rows(0), num_cols(0), nnz(0) {} // 默认
    SparseMatrix(int rows, int cols);
    
    // 从三重格式构造
    SparseMatrix(int rows, int cols, 
                const std::vector<int>& row_indices,
                const std::vector<int>& col_indices,
                const Vector& vals);            // 使用Vector替代

    // 获取矩阵信息
    int getNumRows() const;
    int getNumCols() const;
    size_t getNNZ() const;

    // 矩阵操作
    double get(int row, int col) const;
    void set(int row, int col, double value);
    void add(int row, int col, double value);
    Vector multiply(const Vector& x) const;      // 使用Vector替代

    // 打印矩阵
    void print() const;
};

#endif // SPARSE_MATRIX_HPP