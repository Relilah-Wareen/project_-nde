#ifndef RK_METHOD_HPP
#define RK_METHOD_HPP

#include "ODESolver.hpp"

class RKMethod : public ODESolver {
public:
    RKMethod(int p, double T, std::function<Vector(const Vector&, double)> f);
    void initialize(const Vector& init_value, int n) override;
    void solve() override;
protected:
    virtual Vector computeStep(const Vector& u_n, double t_n, double dt) = 0;
};

class ClassicalRK4 : public RKMethod {
public:
    ClassicalRK4(int p, double T, std::function<Vector(const Vector&, double)> f);
protected:
    Vector computeStep(const Vector& u_n, double t_n, double dt) override;
};


class ESDIRK64 : public RKMethod {
private:
    static const Vector c;
    static const std::vector<Vector> A;
    static const Vector b;

protected:
    Vector computeStep(const Vector& u_n, double t_n, double dt) override;
public:
    ESDIRK64(int p, double T, std::function<Vector(const Vector&, double)> f)
        : RKMethod(p, T, f) {}
};

class GaussLegendre : public RKMethod {
private:
    int s; // 阶段数
    static const std::vector<Vector> c_gauss;
    static const std::vector<std::vector<Vector>> A_gauss;
    static const std::vector<Vector> b_gauss;

protected:
    Vector computeStep(const Vector& u_n, double t_n, double dt) override;

public:
    GaussLegendre(int s, double T, std::function<Vector(const Vector&, double)> f)
        : RKMethod(2 * s, T, f), s(s) {
        if (s != 2 && s != 3 && s != 4 && s != 5) { 
            throw std::invalid_argument("Gauss-Legendre method only supports s=2,3,4,5.");
        }
    }
};


#endif // RK_METHOD_HPP