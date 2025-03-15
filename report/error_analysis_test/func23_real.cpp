#include <fstream>
#include <iostream>
#include <cmath>

// 真解 u(x, y) = exp(x + y)
double true_solution_exp(double x, double y) {
    return exp(x + y);
}

// 真解 u(x, y) = x + y
double true_solution_linear(double x, double y) {
    return x + y;
}

void saveTrueSolution(const std::string& filename, double (*true_solution)(double, double)) {
    // 网格大小
    int grid_num = 64;
    double h = 1.0 / grid_num; // 网格步长

    // 保存结果到 CSV 文件
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return;
    }

    // 计算并输出真解
    for (int row = grid_num; row >= 0; row--) {
        double y = row * h; // 当前行的 y 坐标
        for (int col = 0; col <= grid_num; col++) {
            double x = col * h; // 当前列的 x 坐标
            outFile << true_solution(x, y);
            if (col < grid_num) outFile << ",";
        }
        outFile << std::endl;
    }

    outFile.close();
    std::cout << "True solution saved to " << filename << std::endl;
}

int main() {
    // 计算并保存 u(x, y) = exp(x + y) 的真解
    saveTrueSolution("../../data/true_solution_exp.csv", true_solution_exp);

    // 计算并保存 u(x, y) = x + y 的真解
    saveTrueSolution("../../data/true_solution_linear.csv", true_solution_linear);

    return 0;
}