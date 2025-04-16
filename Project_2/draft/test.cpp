#include "src/BVPSolver.hpp"
#include "src/Operator.hpp"
#include "src/BoundCon.hpp"
#include "src/Function.hpp"
#include <iostream>
#include <cmath>
#include <chrono>  // 高精度计时库

// 测试函数1D real u = exp(sin(x))
double test_func(double x) {
    return (sin(x) - cos(x) * cos(x)) * exp(sin(x));
}

// 测试函数2D real u = exp(sin(x)+y)
double f2D(double x, double y) {
    return -exp(y+sin(x)) + sin(x)*exp(y+sin(x)) - cos(x)*cos(x)*exp(y+sin(x));
}

double Sq_Dirichlet_bound(double x, double y) {
    return exp(y+sin(x));
}

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


int main() {
    try {
        // 创建测试函数
        // Function f(test_func);
        Function f(f2D);
        Function Dbound_2d(Sq_Dirichlet_bound);
        Function N1bound(Sq_Neumann_bound1);
        Function N2bound(Sq_Neumann_bound2);
        Function N3bound(Sq_Neumann_bound3);
        Function N4bound(Sq_Neumann_bound4);
        // 创建操作符对象
        Injection<1> injection1D;
        LinearInterpolation<1> interpolation1D;
        FullWeighting<1> fullweight1D;
        QuadraticInterpolation<1> QuadInter1D;
        Injection<2> injection2D;
        LinearInterpolation<2> interpolation2D;
        FullWeighting<2> fullweight2D;
        QuadraticInterpolation<2> QuadInter2D;
        // 创建边界条件
        // BoundCon bc(1, BCType::Dirichlet, 1, exp(sin(1)));
        // BoundCon bc(1, BCType::Neumann, 1, cos(1)*exp(sin(1)) );
        // BoundCon bc(1, BCType::Mixed, 1, exp(sin(1)), Mixed1DType::ND);
        BoundCon bc(2,BCType::Dirichlet, Dbound_2d);
        // BoundCon bc(2,BCType::Neumann, N1bound,N2bound,N3bound,N4bound);
        // BoundCon bc(2,BCType::Mixed, Mixed2DType::DDNN, Dbound_2d,Dbound_2d,N3bound,N4bound);
        // bc.set_Irregular();

        // // 创建求解器
        int grid_num = 64; // 网格数
        int nu1 = 3;      // 松弛次数1
        int nu2 = 3;      // 松弛次数2
        double epsilon = 1e-8; // 容差
        
        // BVPSolver<1> solver(fullweight1D, QuadInter1D, grid_num, bc, f, nu1, nu2, epsilon);
        BVPSolver<2> solver(injection2D, interpolation2D, grid_num, bc, f, nu1, nu2, epsilon);
        // 构建线性系统
        solver.BuildLinearSystem();
        
        // 保存矩阵和向量（可选） 仅调试程序时用,在验证速度时不要选!!!!
        // solver.saveMatrixToCSV("matrix_A.csv");
        // solver.saveVectorToCSV("vector_F.csv");
        


        // ===== 测试 V-Cycle 求解时间 =====
        // auto start_vc = std::chrono::high_resolution_clock::now(); // 开始计时
        // solver.SolveByVC();  // 调用 V-Cycle 求解
        // auto end_vc = std::chrono::high_resolution_clock::now();   // 结束计时
        // auto duration_vc = std::chrono::duration_cast<std::chrono::microseconds>(end_vc - start_vc);
        // std::cout << "V-Cycle Time: " << duration_vc.count() << " μs" << std::endl;
        
        // ===== 测试 FMG 求解时间 =====
        auto start_fmg = std::chrono::high_resolution_clock::now(); // 开始计时
        solver.SolveByLU();
                // solver.SolveByRelaxOnly(4.0/5);
        // solver.SolveByFMG();  // 调用 FMG 求解
        auto end_fmg = std::chrono::high_resolution_clock::now();   // 结束计时
        auto duration_fmg = std::chrono::duration_cast<std::chrono::microseconds>(end_fmg - start_fmg);
        std::cout << "FMG Time: " << duration_fmg.count() << " μs" << std::endl;
        
        // 输出结果
        // solver.normalizeNeumannSolution();
        solver.printResultToCSV("solution.csv");
        // std::cout << "Iteration count: " << solver.getCount() << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}