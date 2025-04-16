#ifndef BVP_SOLVER_HPP
#define BVP_SOLVER_HPP

#include <vector>
#include <functional>
#include <array>
#include <variant>
#include <iostream>
#include "SparseMatrix.hpp"

enum class BCType { Dirichlet, Neumann };

// 一维边界组合 (左,右)
enum class Mixed1DType { DD, DN, ND, NN };

// 二维边界组合 (左,右,下,上)
enum class Mixed2DType {
    DDDD, DDDN, DDND, DDNN,
    DNDD, DNDN, DNND, DNNN,
    NDDD, NDDN, NDND, NDNN,
    NNDD, NNDN, NNND, NNNN
};

template<int Dim>
class BVPSolver {
public:
    // 1D 构造函数
    BVPSolver(int n, Mixed1DType bc_type, double left_val, double right_val);
    
    // 2D 构造函数
    BVPSolver(int n, BCType bc_type, std::function<double(double, double)> bc_func);
    BVPSolver(int n, Mixed2DType bc_type, std::array<std::function<double(double, double)>, 4> bc_funcs);

    // 设置右端项
    void setRightHandSide(std::function<double(double)> f); // 1D
    void setRightHandSide(std::function<double(double, double)> f); // 2D

    void buildLinearSystem();
    void printSystem() const;
    std::vector<double> solve();

private:
    int n;
    double h;
    SparseMatrix A;
    std::vector<double> rhs;

    // 右端项存储
    std::variant<
        std::function<double(double)>,
        std::function<double(double, double)>
    > f;

    // 边界条件存储
    struct Boundary1D { /* ... */ };
    struct Boundary2D { /* ... */ };
    std::variant<Boundary1D, Boundary2D> boundary;

    void initialize();
};

#endif // BVP_SOLVER_HPP