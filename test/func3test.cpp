#include "../src/BVPSolver.hpp"
#include "../src/ProblemDomain.hpp"
#include <fstream>
#include <iostream>

// 源函数 f(x, y) = -exp(x + y)
double f(double x, double y) {
    return -exp(x + y);
}

// Dirichlet 边界条件：u(x, y) = exp(x + y)
double Sq_Dirichlet_bound(double x, double y) {
    return exp(x + y);
}

// 圆形边界条件：u(x, y) = exp(x + y)
double Circle_Dirichlet_bound(double x, double y) {
    return exp(x + y);
}

int main() {
    // 网格大小
    int grid_num = 64;

    // 创建 ProblemDomain，设置圆形区域
    double centerX = 0.5; // 圆心 x 坐标
    double centerY = 0.5; // 圆心 y 坐标
    double radius = 0.2;  // 圆半径
    ProblemDomain domain(DomainType::Circular, grid_num, centerX, centerY, radius);

    // 设置方形边界条件（Dirichlet）
    std::array<std::function<double(double, double)>, 4> boundary_funcs = {
        Sq_Dirichlet_bound,  // 左边界
        Sq_Dirichlet_bound,  // 右边界
        Sq_Dirichlet_bound,  // 下边界
        Sq_Dirichlet_bound   // 上边界
    };
    domain.setSquareBoundaryCondition(BCType::Dirichlet, boundary_funcs);

    // 设置圆形边界条件（Dirichlet）
    domain.setCircleBoundaryCondition(BCType::Dirichlet, Circle_Dirichlet_bound);

    // 设置 (0,0) 处的值（对于 Dirichlet 边界条件，这个值不会被使用）
    domain.SetValueAt00(exp(0.0 + 0.0)); // u(0,0) = exp(0 + 0) = 1

    // 创建 BVPSolver 并求解
    BVPSolver solver(domain, f);
    solver.BuildLineraSystem();
    std::vector<double> X = solver.solve();

    // 保存结果到 CSV 文件
    std::ofstream outFile("../data/func3.csv");
    if (!outFile) {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return 1;
    }

    // 输出结果到 CSV 文件
    for (int row = grid_num; row >= 0; row--) {
        for (int col = 0; col < grid_num + 1; col++) {
            outFile << X[row * (grid_num + 1) + col];
            if (col < grid_num) outFile << ",";
        }
        outFile << std::endl;
    }

    outFile.close();
    std::cout << "Results saved to func3.csv" << std::endl;

    return 0;
}