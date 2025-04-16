#ifndef BVPSOLVER_HPP
#define BVPSOLVER_HPP

#include "SparseMatrix.hpp"
#include "Function.hpp"
#include "BoundCon.hpp"
#include "Operator.hpp"
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
typedef SparseMatrix Matrix;

template <int Dim>
class BVPSolver {
private:
    const Operator<Dim>& restriction;
    const Operator<Dim>& prolongation;
    std::map<int, Matrix> A;
    Vector F;
    Vector result;
    BoundCon boundary_cond;  // 直接存储边界条件对象
    Function func;
    int GridNum;
    int nu1;
    int nu2;
    double err;
    int cnt; // 内置计数器
    // Vector restriction_op(const int& n, const Vector& r);
    // Vector interpolation(const int& n, const Vector& r);
    Vector VC(const int& n, Vector v, const Vector& f);
    Vector FMG(const int& n, const Vector& f);
    Vector relax(const Matrix& A, const Vector& v, const Vector& f, double omega);
    Vector residual(const Matrix& A, const Vector& v, const Vector& f) const;
    Vector getZeroVector(int n) const;
    Vector LUSolver(const std::vector<std::vector<double>>& A, const Vector& b) const;  
public:
    BVPSolver(const Operator<Dim>& restrict, 
              const Operator<Dim>& prolong,
              const int& n,
              const BoundCon& bc,  // 直接传入BoundCon对象
              const Function& f,
              const int& nu1,
              const int& nu2,
              const double& epsilon);
    Vector plus_vec(const Vector& v1, const Vector& v2);
    void BuildLinearSystem();
    void BuildForIrr();
    const Vector& getResult() const;
    void printA() const;
    void printF() const;
    void SolveByLU();
    void SolveByRelaxOnly(double omega);
    void SolveByVC();
    void SolveByFMG();
    void printResult() const;
    int getCount() const;
    void saveMatrixToCSV(const std::string& filename) const;
    void saveVectorToCSV(const std::string& filename) const;
    int idx_2d(int n, int i, int j) const { return i + (n + 1) * j; }
    void printResultToCSV(const std::string& filename) const;
    void normalizeNeumannSolution();
};

#endif // BVPSOLVER_HPP