#ifndef BVPSOLVER_HPP
#define BVPSOLVER_HPP

#include "ProblemDomain.hpp"
#include <lapacke.h>
#include <cmath>
#include <vector>
#include <memory>
#include <fstream>

class BVPSolver {
private:
    ProblemDomain pd;
    int GridNum;
    double h;
    std::vector<std::vector<double>> A;
    std::vector<double> F;
    double (*function)(double, double);

public:
    BVPSolver(ProblemDomain pd, double (*function)(double, double));

    void ShowMatrix();
    void ShowVectorF();
    void ShowMatrix(const std::string& filename);
    void ShowVectorF(const std::string& filename);

    void BuildLineraSystem();
    std::vector<double> solve();
    int getGridNum() const;
};

#endif // BVPSOLVER_HPP