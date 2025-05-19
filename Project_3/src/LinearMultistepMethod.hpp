#ifndef LINEAR_MULTISTEP_METHOD_HPP
#define LINEAR_MULTISTEP_METHOD_HPP

#include<iostream>
#include <vector>
#include <functional>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <fstream>
#include "ODESolver.hpp"


class LinearMultistepMethod : public ODESolver {
protected:
    int s; // 多步法的步数s
public:
    virtual ~LinearMultistepMethod() = default;
    LinearMultistepMethod(int p, double T, std::function<Vector(const Vector&, double)> f);
    Vector RK_4onestep(const Vector& u_n, double t); // 用于初值的单步RK4
    void initialize(const Vector& init_value, int n) override;
    virtual void solve() = 0;
};
    
class AdamsMethod : public LinearMultistepMethod {
protected:
    virtual double beta(int i, int j) const = 0;
    AdamsMethod(const int& p, const double& T, const std::function<Vector(const Vector&, double)>& f);
};

class AdamsBashforth : public AdamsMethod {
private:
    static const std::unordered_map<int, Vector> beta_coeffs;
    double beta(int i, int j) const override; 
protected:
    using LinearMultistepMethod::s;
    using LinearMultistepMethod::Dim;
public:
    AdamsBashforth(const int& p, const double& T, const std::function<Vector(const Vector&, double)>& f);
    void solve() override;
};

class AdamsMoulton : public AdamsMethod {
private:
    static const std::unordered_map<int, Vector> beta_coeffs;
    double beta(int i, int j) const override; 
protected:
    using LinearMultistepMethod::s;
    using LinearMultistepMethod::Dim;
public:
    AdamsMoulton(const int& p, const double& T, const std::function<Vector(const Vector&, double)>& f);
    void solve() override;
};

class BDFMethod : public LinearMultistepMethod {
private:
    static const std::unordered_map<int, Vector> alpha_coeffs;
    static const std::unordered_map<int, double> beta_s;
    double alpha(int i, int j);
    double beta(int s);
public:
    BDFMethod(const int& p, const double& T, const std::function<Vector(const Vector&, double)>& f);
    void solve() override;
};

#endif // LINEAR_MULTISTEP_METHOD_HPP