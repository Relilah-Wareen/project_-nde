#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include <vector>
#include <functional>
#include <string>
#include <stdexcept>
#include <fstream>
#include <cmath>

typedef std::vector<double> Vector;

// 重载数乘运算符：double * Vector
inline Vector operator*(double scalar, const Vector& vec) {
    Vector result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = scalar * vec[i];
    }
    return result;
}

// 重载数乘运算符：Vector * double
inline Vector operator*(const Vector& vec, double scalar) {
    return scalar * vec; // 复用上面的实现
}

inline Vector operator+(const Vector& v1, const Vector& v2) {
    if (v1.size() != v2.size()) throw std::invalid_argument("Vector sizes mismatch");
    Vector result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

inline Vector operator-(const Vector& v1, const Vector& v2) {
    if (v1.size() != v2.size()) throw std::invalid_argument("Vector sizes mismatch");
    Vector result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}


inline Vector& operator+=(Vector& v1, const Vector& v2) {
    if (v1.size() != v2.size()) throw std::invalid_argument("Vector sizes mismatch");
    for (size_t i = 0; i < v1.size(); ++i) {
        v1[i] += v2[i];
    }
    return v1;
}

inline double Inf_norm(const Vector& vec) {
    if (vec.empty()) {
        throw std::invalid_argument("Vector is empty.");
    }

    return *std::max_element(vec.begin(), vec.end(), 
                             [](double a, double b) { return std::abs(a) < std::abs(b); });
}

inline double L2_norm(const Vector& vec) {
    if (vec.empty()) {
        throw std::invalid_argument("Vector is empty.");
    }

    double sum_of_squares = 0.0;
    for (const auto& elem : vec) {
        sum_of_squares += elem * elem;
    }

    return std::sqrt(sum_of_squares);
}


class ODESolver {
protected:
    int p;          // 阶数
    double T;       // 总时间
    int Num;        // 总步数
    double k;       // 步长
    int Dim;        // 解向量的维度
    std::function<Vector(const Vector&, double)> F; // 右端项函数
    std::vector<Vector> result; // 数值解存储

public:
    virtual ~ODESolver() = default;
    ODESolver(int p, double T, std::function<Vector(const Vector&, double)> f);
    virtual void initialize(const Vector& init_value, int n) = 0;
    virtual void solve() = 0;
    const std::vector<Vector>& get_results() const;
    void printforU1U2(const std::string& filename = "u1u2.csv") const;
    double computePeriodError();
};

#endif // ODE_SOLVER_HPP