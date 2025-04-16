#include "Operator.hpp"
#include <iostream>
#include <iomanip>

void print_vector(const std::vector<double>& v, const std::string& name) {
    std::cout << name << ": [";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << std::fixed << std::setprecision(2) << v[i];
        if (i != v.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}

void test_injection() {
    Injection<1> inj;
    
    std::vector<double> v1 = {1.0, 2.0, 3.0, 4.0, 5.0,7,8,9};
    auto result1 = inj(v1);
    print_vector(v1, "Input ");
    print_vector(result1, "Output");
    std::cout << std::endl;
    
    std::vector<double> v2 = {10.0, 20.0, 30.0, 50, 60};
    auto result2 = inj(v2);
    print_vector(v2, "Input ");
    print_vector(result2, "Output");
    std::cout << std::endl;
}

void test_interpolation() {
    LinearInterpolation<1> interp;
    
    std::vector<double> v1 = {1.0, 3.0, 5.0};
    auto result1 = interp(v1);
    print_vector(v1, "Input ");
    print_vector(result1, "Output");
    std::cout << std::endl;
    
    std::vector<double> v2 = {10.0, 20.0, 40,50,60};
    auto result2 = interp(v2);
    print_vector(v2, "Input ");
    print_vector(result2, "Output");
    std::cout << std::endl;
}

int main() {
    std::cout << "=== Testing Injection ===" << std::endl;
    test_injection();
    
    std::cout << "=== Testing Linear Interpolation ===" << std::endl;
    test_interpolation();
    
    return 0;
}