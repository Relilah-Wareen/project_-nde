#include "BVPSolver.hpp"
#include <iostream>

void test_1d_system() {
    const int n = 8; 
    
    // 创建求解器（两边Dirichlet）
    BVPSolver<1> solver(n, Mixed1DType::DD, 1.0, 2.0);
    
    // 设置右端项 f(x) = 10
    solver.setRightHandSide([](double x) { return 10.0; });
    
    // 构建并输出系统
    solver.buildLinearSystem();
    solver.printSystem();
}

int main() {
    test_1d_system();
    return 0;
}