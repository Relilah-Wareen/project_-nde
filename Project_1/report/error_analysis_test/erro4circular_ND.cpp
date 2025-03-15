#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

// 读取 CSV 文件中的数据
std::vector<std::vector<int>> readBoundaryPoints(const std::string& filename, int num) {
    std::vector<std::vector<int>> boundaryPoints(num + 1, std::vector<int>(num + 1, 0));
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Error: Could not open file " << filename << " for reading." << std::endl;
        return boundaryPoints;
    }

    for (int j = 0; j <= num; j++) {
        for (int i = 0; i <= num; i++) {
            int value;
            inFile >> value;
            boundaryPoints[j][i] = value;
            if (inFile.peek() == ',') inFile.ignore();
        }
    }

    inFile.close();
    return boundaryPoints;
}

// 读取数值解或精确解文件
std::vector<std::vector<double>> readSolution(const std::string& filename, int num) {
    std::vector<std::vector<double>> solution(num + 1, std::vector<double>(num + 1, 0.0));
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Error: Could not open file " << filename << " for reading." << std::endl;
        return solution;
    }

    for (int j = 0; j <= num; j++) {
        for (int i = 0; i <= num; i++) {
            double value;
            inFile >> value;
            solution[j][i] = value;
            if (inFile.peek() == ',') inFile.ignore();
        }
    }

    inFile.close();
    return solution;
}

// 计算边界点的 L2 误差
double calculateBoundaryL2Error(const std::vector<std::vector<int>>& boundaryPoints,
                                const std::vector<std::vector<double>>& numerical,
                                const std::vector<std::vector<double>>& exact,
                                int num) {
    double h = 1.0 / num;
    double l2_error = 0.0;

    // 遍历所有格点
    for (int j = 0; j <= num; j++) {
        for (int i = 0; i <= num; i++) {
            // 如果是边界点
            if (boundaryPoints[j][i] == 1) {
                double error = exact[j][i] - numerical[j][i];
                l2_error += error * error;
            }
        }
    }

    l2_error = std::sqrt(l2_error * h * h);
    return l2_error;
}

int main() {
    // 网格大小列表
    std::vector<int> grid_sizes = {8, 16, 32, 64};

    // 遍历每个网格大小
    for (int num : grid_sizes) {
        // 读取边界点文件
        std::string boundary_file = "../../data/boundary_points_N" + std::to_string(num) + ".csv";
        auto boundaryPoints = readBoundaryPoints(boundary_file, num);

        // 读取数值解文件
        std::string numerical_file = "../../data/result_Circular_Neumann_Dirichlet_N" + std::to_string(num) + ".csv";
        auto numerical = readSolution(numerical_file, num);

        // 读取精确解文件
        std::string exact_file = "../../data/realnum_" + std::to_string(num) + ".csv";
        auto exact = readSolution(exact_file, num);

        // 计算边界点的 L2 误差
        double l2_error = calculateBoundaryL2Error(boundaryPoints, numerical, exact, num);

        // 输出结果
        std::cout << "Grid size: " << num << ", Boundary L2 Error: " << std::setprecision(10) << l2_error << std::endl;
    }

    return 0;
}