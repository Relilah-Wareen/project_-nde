#include "BVPSolver.hpp"
#include <stdexcept>
#include <iostream>

template<int Dim>
BVPSolver<Dim>::BVPSolver(int n, Mixed1DType bc_type, double left_val, double right_val)
    : n(n), h(1.0 / n), A(n + 1, n + 1), boundary{Boundary1D{bc_type, left_val, right_val}} 
{
    if (Dim != 1) {
        throw std::runtime_error("This constructor is for 1D only");
    }
    initialize();
}

template<int Dim>
BVPSolver<Dim>::BVPSolver(int n, BCType bc_type, std::function<double(double, double)> bc_func)
    : n(n), h(1.0 / n), A((n + 1) * (n + 1), (n + 1) * (n + 1)), boundary{Boundary2D{bc_type, bc_func}} 
{
    if (Dim != 2) {
        throw std::runtime_error("This constructor is for 2D only");
    }
    initialize();
}

template<int Dim>
void BVPSolver<Dim>::setRightHandSide(std::function<double(double)> f) {
    if (Dim != 1) {
        throw std::runtime_error("This method is for 1D only");
    }
    this->f = f;
}

template<int Dim>
void BVPSolver<Dim>::buildLinearSystem() {
    if (Dim == 1) {
        // 1D 离散化代码...
    } 
    else if (Dim == 2) {
        // 2D 离散化代码...
    }
    else {
        throw std::runtime_error("Invalid dimension");
    }
}

// 其他方法类似修改...