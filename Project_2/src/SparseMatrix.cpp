#include "SparseMatrix.hpp"
#include <stdexcept>
#include <iostream>

SparseMatrix::SparseMatrix(int rows, int cols) 
    : num_rows(rows), num_cols(cols), nnz(0) {
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }
    row_ptr.resize(rows + 1, 0);
}

SparseMatrix::SparseMatrix(int rows, int cols, 
                         const std::vector<int>& row_indices,
                         const std::vector<int>& col_indices,
                         const Vector& vals) 
    : num_rows(rows), num_cols(cols) {
    if (row_indices.size() != col_indices.size() || row_indices.size() != vals.size()) {
        throw std::invalid_argument("Input vectors must have the same size");
    }
    
    // 统计每行非零元素个数
    row_ptr.resize(num_rows + 1, 0);
    for (size_t i = 0; i < row_indices.size(); ++i) {
        if (row_indices[i] < 0 || row_indices[i] >= num_rows ||
            col_indices[i] < 0 || col_indices[i] >= num_cols) {
            throw std::out_of_range("Index out of matrix bounds");
        }
        row_ptr[row_indices[i] + 1]++;
    }
    
    // 计算行指针
    for (int i = 1; i <= num_rows; ++i) {
        row_ptr[i] += row_ptr[i - 1];
    }
    
    // 填充values和columns
    nnz = row_ptr[num_rows];
    values.resize(nnz);
    columns.resize(nnz);
    
    std::vector<int> row_counter(num_rows, 0);
    for (size_t i = 0; i < row_indices.size(); ++i) {
        int row = row_indices[i];
        int pos = row_ptr[row] + row_counter[row];
        values[pos] = vals[i];
        columns[pos] = col_indices[i];
        row_counter[row]++;
    }
}

int SparseMatrix::getNumRows() const { return num_rows; }
int SparseMatrix::getNumCols() const { return num_cols; }
size_t SparseMatrix::getNNZ() const { return nnz; }

double SparseMatrix::get(int row, int col) const {
    if (row < 0 || row >= num_rows || col < 0 || col >= num_cols) {
        throw std::out_of_range("Index out of bounds");
    }
    
    for (int i = row_ptr[row]; i < row_ptr[row + 1]; ++i) {
        if (columns[i] == col) {
            return values[i];
        }
    }
    return 0.0;
}

void SparseMatrix::set(int row, int col, double value) {
    if (row < 0 || row >= num_rows || col < 0 || col >= num_cols) {
        throw std::out_of_range("Index out of bounds");
    }
    
    for (int i = row_ptr[row]; i < row_ptr[row + 1]; ++i) {
        if (columns[i] == col) {
            values[i] = value;
            return;
        }
    }
    
    if (value != 0.0) {
        insert(row, col, value);
    }
}

void SparseMatrix::add(int row, int col, double value) {
    if (row < 0 || row >= num_rows || col < 0 || col >= num_cols) {
        throw std::out_of_range("Index out of bounds");
    }
    
    // 查找元素是否存在
    for (int i = row_ptr[row]; i < row_ptr[row + 1]; ++i) {
        if (columns[i] == col) {
            values[i] += value;
            return;
        }
    }
    
    // 如果元素不存在且value不为0，则插入新元素
    if (value != 0.0) {
        insert(row, col, value);
    }
}

Vector SparseMatrix::multiply(const Vector& x) const {
    if (x.size() != static_cast<size_t>(num_cols)) {
        throw std::invalid_argument("Vector size does not match matrix columns");
    }
    
    Vector y(num_rows, 0.0);
    for (int row = 0; row < num_rows; ++row) {
        for (int i = row_ptr[row]; i < row_ptr[row + 1]; ++i) {
            y[row] += values[i] * x[columns[i]];
        }
    }
    return y;
}

void SparseMatrix::print() const {
    for (int row = 0; row < num_rows; ++row) {
        for (int col = 0; col < num_cols; ++col) {
            std::cout << get(row, col) << "\t";
        }
        std::cout << std::endl;
    }
}

void SparseMatrix::insert(int row, int col, double value) {
    int insert_pos = row_ptr[row + 1];
    values.insert(values.begin() + insert_pos, value);
    columns.insert(columns.begin() + insert_pos, col);
    for (int i = row + 1; i <= num_rows; ++i) {
        row_ptr[i]++;
    }
    nnz++;
}