#include "../src/BVPSolver.hpp"
#include "../src/ProblemDomain.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <map>

double f(double x, double y) {
    return -exp(y+sin(x)) + sin(x)*exp(y+sin(x)) - cos(x)*cos(x)*exp(y+sin(x));
}
 
double Sq_Dirichlet_bound(double x, double y) {
    return exp(y+sin(x));
}

// 顺序依次为左右下上
double Sq_Neumann_bound1(double x, double y){
    return -exp(y);
}

double Sq_Neumann_bound2(double x, double y){
    return exp(y+sin(1))*cos(1);
}

double Sq_Neumann_bound3(double x, double y){
    return -exp(sin(x));
}

double Sq_Neumann_bound4(double x, double y){
    return exp(1+sin(x));
}


double ftest(double x, double y) {
    return x+y;
}

double Sq_Neumann_x(double x, double y){
    return cos(x)*exp(y+sin(x));
}

double Sq_Neumann_y(double x, double y){
    return exp(y+sin(x));
}


using json = nlohmann::json;

// 定义函数映射表
std::map<std::string, std::function<double(double, double)>> function_map = {
    {"Sq_Dirichlet_bound", Sq_Dirichlet_bound},
    {"Sq_Neumann_bound1", Sq_Neumann_bound1},
    {"Sq_Neumann_bound2", Sq_Neumann_bound2},
    {"Sq_Neumann_bound3", Sq_Neumann_bound3},
    {"Sq_Neumann_bound4", Sq_Neumann_bound4}
};

int main() {
    // 读取 JSON 文件
    std::ifstream config_file("../json/dirichlet_test.json");
    if (!config_file.is_open()) {
        std::cerr << "Error: Could not open config file." << std::endl;
        return 1;
    }

    json config;
    config_file >> config;

    // 遍历测试用例
    for (const auto& test_case : config["test_cases"]) {
        std::string name = test_case["name"];
        int grid_num = test_case["grid_num"];
        double value_at_00 = test_case["value_at_00"];

        // 创建 ProblemDomain
        ProblemDomain domain(DomainType::Default, grid_num);

        // 设置方形边界条件
        std::array<std::function<double(double, double)>, 4> boundary_funcs;
        for (int i = 0; i < 4; ++i) {
            std::string func_name = test_case["square_boundary"]["functions"][i];
            boundary_funcs[i] = function_map[func_name];
        }
        domain.setSquareBoundaryCondition(BCType::Dirichlet, boundary_funcs);

        // 设置 (0,0) 处的值
        domain.SetValueAt00(value_at_00);

        // 创建 BVPSolver 并求解
        BVPSolver solver(domain, f);
        solver.BuildLineraSystem();
        std::vector<double> X = solver.solve();

        // 保存结果
        std::ofstream outFile("../data/result_" + name + ".csv");
        if (!outFile) {
            std::cerr << "Error: Could not open file for writing." << std::endl;
            return 1;
        }

        for (int row = grid_num; row >= 0; row--) {
            for (int col = 0; col < grid_num + 1; col++) {
                outFile << X[row * (grid_num + 1) + col] << ",";
            }
            outFile << std::endl;
        }

        outFile.close();
        std::cout << "Test case " << name << " completed. Results saved to result_" << name << ".csv" << std::endl;
    }

    return 0;
}