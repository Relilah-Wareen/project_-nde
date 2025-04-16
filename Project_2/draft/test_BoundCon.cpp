#include "BoundCon.hpp"
#include <iostream>
#include <cmath>

// 辅助函数：打印测试结果
void print_test_result(const std::string& name, bool passed) {
    std::cout << (passed ? "[PASS] " : "[FAIL] ") << name << std::endl;
}

// 测试1D边界条件
void test_1d_boundaries() {
    try {
        // 测试Dirichlet边界
        BoundCon bc_dirichlet(1, BCType::Dirichlet, 1.0, 2.0);
        bool dirichlet_ok = (bc_dirichlet.getLeft1D() == 1.0) && 
                           (bc_dirichlet.getRight1D() == 2.0);
        print_test_result("1D Dirichlet Boundary", dirichlet_ok);

        // 测试Mixed边界 (DN类型)
        BoundCon bc_mixed(1, BCType::Mixed, 3.0, 4.0, Mixed1DType::DN);
        bool mixed_ok = (bc_mixed.getLeft1D() == 3.0) && 
                       (bc_mixed.getRight1D() == 4.0);
        print_test_result("1D Mixed Boundary (DN)", mixed_ok);

        // 测试错误处理（错误维度）
        try {
            BoundCon bc_invalid(3, BCType::Dirichlet, 0.0, 0.0);
            print_test_result("1D Invalid Dimension", false);
        } catch (const std::invalid_argument&) {
            print_test_result("1D Invalid Dimension", true);
        }
    } catch (...) {
        print_test_result("1D Tests", false);
    }
}

// 测试2D边界条件
void test_2d_boundaries() {
    // 定义测试函数
    Function zero_func([](double x, double y) { return 0.0; });
    Function left_func([](double x, double y) { return y; });      // u(0,y) = y
    Function right_func([](double x, double y) { return 1.0-y; }); // u(1,y) = 1-y
    Function bottom_func([](double x, double y) { return x; });    // u(x,0) = x
    Function top_func([](double x, double y) { return 1.0-x; });   // u(x,1) = 1-x

    try {
        // 测试Uniform Dirichlet边界
        BoundCon bc_uniform(2, BCType::Dirichlet, zero_func);
        bool uniform_ok = (bc_uniform.getBound2D(0.0, 0.5) == 0.0) &&
                          (bc_uniform.getBound2D(0.3, 1.0) == 0.0);
        print_test_result("2D Uniform Dirichlet", uniform_ok);

        // 测试Mixed边界 (DNDN类型)
        BoundCon bc_mixed(2, BCType::Mixed, Mixed2DType::DNDN,
                         left_func, right_func, bottom_func, top_func);
        bool mixed_ok = (bc_mixed.getBound2D(0.0, 0.5) == 0.0) &&    // 左边
                        (bc_mixed.getBound2D(1.0, 0.2) == 0.8) &&    // 右边
                        (bc_mixed.getBound2D(0.3, 0.0) == 0.3) &&    // 下边
                        (bc_mixed.getBound2D(0.7, 1.0) == 0.3);      // 上边
        print_test_result("2D Mixed Boundary (DNDN)", mixed_ok);

        // 测试角点
        bool corner_ok = (bc_mixed.getBound2D(0.0, 0.0) == 0.0) &&   // 左下角
                         (bc_mixed.getBound2D(1.0, 1.0) == 0.0);     // 右上角
        print_test_result("2D Corner Points", corner_ok);

        // 测试错误处理（内部点）
        try {
            bc_mixed.getBound2D(0.5, 0.5);
            print_test_result("2D Internal Point Check", false);
        } catch (const std::runtime_error&) {
            print_test_result("2D Internal Point Check", true);
        }
    } catch (...) {
        print_test_result("2D Tests", false);
    }
}

int main() {
    std::cout << "=== Testing 1D Boundaries ===" << std::endl;
    test_1d_boundaries();

    std::cout << "\n=== Testing 2D Boundaries ===" << std::endl;
    test_2d_boundaries();

    return 0;
}