#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

// 读取 CSV 文件中的数据
std::vector<double> readCSV(const std::string& filename, int num) {
    std::vector<double> data((num + 1) * (num + 1));
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Error: Could not open file " << filename << " for reading." << std::endl;
        return data;
    }

    for (int row = num; row >= 0; row--) {
        for (int col = 0; col <= num; col++) {
            double value;
            inFile >> value;
            data[row * (num + 1) + col] = value;
            if (inFile.peek() == ',') inFile.ignore();
        }
    }

    inFile.close();
    return data;
}

// 计算内部格点的 L2 误差
double calculateL2Error(const std::vector<double>& numerical, const std::vector<double>& exact, int num) {
    double h = 1.0 / num;
    double l2_error = 0.0;

    // 遍历内部格点（跳过边界）
    for (int j = 1; j < num; j++) {
        for (int i = 1; i < num; i++) {
            double error = exact[j * (num + 1) + i] - numerical[j * (num + 1) + i];
            l2_error += error * error;
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
        // 读取数值解和精确解
        std::string numerical_file = "../data/result_Dirichlet_N" + std::to_string(num) + ".csv";
        std::string exact_file = "../data/realnum_" + std::to_string(num) + ".csv";

        std::vector<double> numerical = readCSV(numerical_file, num);
        std::vector<double> exact = readCSV(exact_file, num);

        // 计算 L2 误差
        double l2_error = calculateL2Error(numerical, exact, num);
        std::cout << "Grid size: " << num << ", L2 Error: " << std::setprecision(10) << l2_error << std::endl;
    }

    return 0;
}