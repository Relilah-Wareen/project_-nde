#define _USE_MATH_DEFINES
#include "../src/IVPFactory.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>
#include <nlohmann/json.hpp>

namespace fs = std::filesystem;
using json = nlohmann::json;

// 确保日志和数据目录存在
void check_directories() {
    fs::create_directories("../data");
    fs::create_directories("../log");
}

Vector rhs_orbit(const Vector& u, double t) {
    const double mu = 0.012277471;
    const double r1_sq = pow(u[0] + mu - 1, 2) + pow(u[1], 2) + pow(u[2], 2);
    const double r2_sq = pow(u[0] + mu, 2) + pow(u[1], 2) + pow(u[2], 2);
    const double r1_cubed = pow(r1_sq, 1.5);
    const double r2_cubed = pow(r2_sq, 1.5);

    return {
        u[3],
        u[4],
        u[5],
        2 * u[4] + u[0] - (mu * (u[0] + mu - 1)) / r1_cubed - ((1 - mu) * (u[0] + mu)) / r2_cubed,
        -2 * u[3] + u[1] - (mu * u[1]) / r1_cubed - ((1 - mu) * u[1]) / r2_cubed,
        -(mu * u[2]) / r1_cubed - ((1 - mu) * u[2]) / r2_cubed
    };
}

void run_test_for_model(const json& model_config) {
    auto& factory = IVPFactory::getInstance();
    std::ofstream logfile("../log/runlog.txt", std::ios::app);

    std::string model_name = model_config["name"];
    Vector u0 = model_config["u0"].get<Vector>();
    double T = model_config["T"];

    for (const auto& method_cfg : model_config["method_configs"]) {
        std::string method = method_cfg["method"];
        const auto& p_n_map = method_cfg["p_n"];

        for (const auto& [p_str, n_value] : p_n_map.items()) {
            try {
                auto start = std::chrono::high_resolution_clock::now();

                int p = std::stoi(p_str);
                int n = n_value;

                // 特殊处理自适应方法（n=0时忽略）
                if (method == "Fehlberg45" || method == "DormandPrince54") {
                    n = 0; // 自适应方法不使用固定步数
                }

                auto solver = factory.createSolver(method, p, T, rhs_orbit);
                solver->initialize(u0, n);
                solver->solve();

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;

                // 记录日志
                logfile << model_name << "," << method << ",p" << p 
                       << ",n" << n << "," << elapsed.count() << "s\n";

                // 生成文件名
                std::string filename = "../data/" + model_name + "_" + method 
                    + "_p" + std::to_string(p) 
                    + "_n" + std::to_string(n) + ".csv";
                solver->printforU1U2(filename);
                std::cout << "Generated: " << filename << std::endl;

            } catch (const std::exception& e) {
                std::cerr << "[ERROR] " << model_name << " " << method 
                        << " p=" << p_str << ": " << e.what() << std::endl;
            }
        }
    }
}

void test_from_json(const std::string& json_path) {
    std::ifstream f(json_path);
    json config = json::parse(f);

    for (const auto& model : config["models"]) {
        run_test_for_model(model);
    }
}

int main() {
    check_directories();
    std::ofstream("../log/runlog.txt", std::ios::trunc); // 清空旧日志
    RegisterAllSolvers();
    test_from_json("../json/models.json");
    return 0;
}
