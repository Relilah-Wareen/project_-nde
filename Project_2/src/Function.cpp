#include "Function.hpp"

Function::Function(std::function<double(double)> f1d)
    : dimension(1), func(f1d) {}

Function::Function(std::function<double(double, double)> f2d)
    : dimension(2), func(f2d) {}

double Function::operator()(double x) const {
    validate_dimension(1);
    return std::get<std::function<double(double)>>(func)(x);
}

double Function::operator()(double x, double y) const {
    validate_dimension(2);
    return std::get<std::function<double(double, double)>>(func)(x, y);
}

void Function::validate_dimension(int expected) const {
    if (dimension != expected) {
        throw std::runtime_error(
            "Function dimension mismatch. Expected " + 
            std::to_string(expected) + "D but has " + 
            std::to_string(dimension) + "D"
        );
    }
}