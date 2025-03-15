#include "../../src/ProblemDomain.hpp"
#include <iostream>

int main() {
    // 网格大小列表
    std::vector<int> grid_sizes = {8, 16, 32, 64};

    // 遍历每个网格大小
    for (int num : grid_sizes) {
        // 创建 ProblemDomain 对象
        ProblemDomain domain(DomainType::Circular, num, 0.5, 0.5, 0.2);

        // 导出边界点到 CSV 文件
        std::string filename = "../../data/boundary_points_N" + std::to_string(num) + ".csv";
        domain.exportBoundaryPointsToCSV(filename);

        std::cout << "Exported boundary points for N = " << num << " to " << filename << std::endl;
    }

    return 0;
}