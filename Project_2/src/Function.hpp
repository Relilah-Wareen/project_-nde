#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <functional>
#include <stdexcept>
#include <variant>

class Function {
public:
    // 构造函数
    Function() = default;
    Function(std::function<double(double)> f1d);
    Function(std::function<double(double, double)> f2d);

    // 调用操作符
    double operator()(double x) const;
    double operator()(double x, double y) const;

    // 获取维度
    int getDimension() const { return dimension; }

    // 检查是否可调用
    bool is1D() const { return dimension == 1; }
    bool is2D() const { return dimension == 2; }

private:
    int dimension;
    std::variant<std::function<double(double)>, 
                 std::function<double(double, double)>> func;

    void validate_dimension(int expected) const;
};

#endif // FUNCTION_HPP