// OperatorTest.cpp
#include "Operator.hpp"
#include <iostream>
#include <vector>

// 打印向量内容
void print_vector(const std::vector<double>& v, const std::string& name) {
    std::cout << name << ": [";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << v[i];
        if (i != v.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}

void test_operators() {
    std::cout << "=== Testing Operators with size=7 and 15 ===" << std::endl;
    
    Injection<1> inject;
    LinearInterpolation<1> interpolate;

    // 测试 size=7 (k=3)
    {
        Vector v7 = {1, 2, 3, 4, 5, 6, 7};  // 简单赋值
        std::cout << "\nCase: size=7 (2^3-1)" << std::endl;
        
        // 测试 Injection
        Vector inj7 = inject(v7);
        print_vector(v7, "Input ");
        print_vector(inj7, "Injection");
        std::cout << "Expected Injection: [2, 4, 6]" << std::endl;
        
        // 测试 LinearInterpolation
        Vector int7 = interpolate(v7);
        print_vector(int7, "Interpolation");
        std::cout << "First 5 of expected Interpolation: [0.5, 1, 1.5, 2, 2.5, ...]" << std::endl;
    }

    // 测试 size=15 (k=4)
    {
        Vector v15 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        std::cout << "\nCase: size=15 (2^4-1)" << std::endl;
        
        // 测试 Injection
        Vector inj15 = inject(v15);
        print_vector(v15, "Input ");
        print_vector(inj15, "Injection");
        std::cout << "Expected Injection: [2, 4, 6, 8, 10, 12, 14]" << std::endl;
        
        // 测试 LinearInterpolation
        Vector int15 = interpolate(v15);
        print_vector(int15, "Interpolation");
        std::cout << "First 5 of expected Interpolation: [0.5, 1, 1.5, 2, 2.5, ...]" << std::endl;
    }
}

int main() {
    test_operators();
    return 0;
}