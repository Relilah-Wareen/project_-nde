#include "Function.hpp"
#include <iostream>
#include <cmath>

// 外部定义的1D函数
double quadratic(double x) {
    return x * x + 2 * x + 1;
}

// 外部定义的2D函数
double trig_product(double x, double y) {
    return std::sin(x) * std::cos(y);
}

int main() {
    // 包装外部定义的1D函数
    std::function<double(double)> f1d_ext = quadratic;
    Function f1(f1d_ext);
    
    // 包装外部定义的2D函数
    std::function<double(double, double)> f2d_ext = trig_product;
    Function f2(f2d_ext);
    
    // 测试1D函数
    std::cout << "=== Testing 1D Function ===" << std::endl;
    std::cout << "f1(2.0) = " << f1(2.0) << std::endl;
    std::cout << "f1(3.0) = " << f1(3.0) << std::endl;
    
    // 测试2D函数
    std::cout << "\n=== Testing 2D Function ===" << std::endl;
    std::cout << "f2(0.0, 0.0) = " << f2(0.0, 0.0) << std::endl;
    std::cout << "f2(M_PI/2, M_PI) = " << f2(M_PI/2, M_PI) << std::endl;
    
    // 测试错误处理
    std::cout << "\n=== Testing Error Handling ===" << std::endl;
    try {
        std::cout << "Attempting to call 1D function with 2 arguments..." << std::endl;
        f1(1.0, 2.0); // 应该抛出异常
    } catch (const std::exception& e) {
        std::cout << "Error caught: " << e.what() << std::endl;
    }
    
    return 0;
}