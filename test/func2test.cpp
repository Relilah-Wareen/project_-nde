#include "../src/BVPSolver.hpp"
#include "../src/ProblemDomain.hpp"
#include <fstream>
#include <iostream>

// 源函数 f(x, y) = 0
double f(double x, double y) {
    return 0.0;
}

// Dirichlet 边界条件：u(x, y) = x + y
double Sq_Dirichlet_bound(double x, double y) {
    return x + y;
}

int main() {
    // 网格大小
    int grid_num = 64;

    // 创建 ProblemDomain
    ProblemDomain domain(DomainType::Default, grid_num);

    // 设置方形边界条件（Dirichlet）
    std::array<std::function<double(double, double)>, 4> boundary_funcs = {
        Sq_Dirichlet_bound,  // 左边界
        Sq_Dirichlet_bound,  // 右边界
        Sq_Dirichlet_bound,  // 下边界
        Sq_Dirichlet_bound   // 上边界
    };
    domain.setSquareBoundaryCondition(BCType::Dirichlet, boundary_funcs);

    // 设置 (0,0) 处的值（对于 Dirichlet 边界条件，这个值不会被使用）
    domain.SetValueAt00(0.0);

    // 创建 BVPSolver 并求解
    BVPSolver solver(domain, f);
    solver.BuildLineraSystem();
    std::vector<double> X = solver.solve();

    // 保存结果到 CSV 文件
    std::ofstream outFile("../data/func2.csv");
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
    std::cout << "Results saved to func2.csv" << std::endl;

    return 0;
}