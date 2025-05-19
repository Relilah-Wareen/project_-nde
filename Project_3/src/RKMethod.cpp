#include "RKMethod.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
RKMethod::RKMethod(int p, double T, std::function<Vector(const Vector&, double)> f)
    : ODESolver(p, T, f) {}

void RKMethod::initialize(const Vector& init_value, int n) {
    Dim = init_value.size();
    Num = n;
    k = T / n;
    result.resize(n + 1, Vector(Dim, 0.0));
    result[0] = init_value;
}

void RKMethod::solve() {
    for (int i = 1; i <= Num; ++i) {
        result[i] = computeStep(result[i-1], (i-1)*k, k);
    }
}

ClassicalRK4::ClassicalRK4(int p, double T, std::function<Vector(const Vector&, double)> f)
    : RKMethod(p, T, f) {
    if (p != 4) {
        throw std::invalid_argument("Classical RK4 requires order 4.");
    }
}



Vector ClassicalRK4::computeStep(const Vector& u_n, double t_n, double dt) {
    Vector k1 = F(u_n, t_n);
    Vector k2 = F(u_n + 0.5 * dt * k1, t_n + 0.5 * dt);
    Vector k3 = F(u_n + 0.5 * dt * k2, t_n + 0.5 * dt);
    Vector k4 = F(u_n + dt * k3, t_n + dt);
    return u_n + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
}

// Butcher tableau coefficients for ESDIRK64
const Vector ESDIRK64::c = {0.0, 0.5, 83.0/250.0, 31.0/50.0, 17.0/20.0, 1.0};

const std::vector<Vector> ESDIRK64::A = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.25, 0.25, 0.0, 0.0, 0.0, 0.0},
    {8611.0/62500.0, -1743.0/31250.0, 0.25, 0.0, 0.0, 0.0},
    {5012029.0/34652500.0, -654441.0/2922500.0, 174375.0/388108.0, 0.25, 0.0, 0.0},
    {1526708209.0/15537626600.0, -71443401.0/120774400.0, 730878875.0/902184768.0, 2285395.0/8070912.0, 0.25, 0.0},
    {82889.0/524892.0, 0.0, 15625.0/83664.0, 69875.0/102672.0, -2260.0/8211.0, 0.25}
};

const Vector ESDIRK64::b = {82889.0/524892.0, 0.0, 15625.0/83664.0, 69875.0/102672.0, -2260.0/8211.0, 0.25};

Vector ESDIRK64::computeStep(const Vector& u_n, double t_n, double dt) {
    int s = c.size();
    std::vector<Vector> y(s);
    std::vector<Vector> k(s);

    double epsilon = 1e-15; // 容差
    int max_iterations = 100; // 最大迭代次数

    for (int i = 0; i < s; ++i) {
        y[i] = u_n;
        for (int j = 0; j < i; ++j) {
            y[i] += dt * A[i][j] * k[j];
        }

        int iterations = 0;
        while (iterations < max_iterations) {
            Vector y_prev = y[i];
            k[i] = F(y_prev, t_n + c[i] * dt);
            y[i] = u_n;
            for (int j = 0; j <= i; ++j) {
                y[i] += dt * A[i][j] * k[j];
            }
            if ((y[i] - y_prev).empty() || std::abs(Inf_norm(y[i] - y_prev)) < epsilon) {
                break;
            }
            iterations++;
        }

        if (iterations >= max_iterations) {
            std::cout << "ESDIRK64: Iteration did not converge within " << max_iterations << " iterations at stage " << i << std::endl;
        }
    }

    Vector y_final = u_n;
    for (int i = 0; i < s; ++i) {
        y_final += dt * b[i] * k[i];
    }

    return y_final;
}

const double w1 = 1.0 / 8 - sqrt(30.0) / 144;
const double w2 = 0.5*sqrt((15 + 2*sqrt(30.0)) / 35);
const double w3 = w2 * (1.0/6 + sqrt(30.0)/24);
const double w4 = w2 * (1.0/21 + 5*sqrt(30.0)/168);
const double w5 = w2 - 2*w3;

const double v1 = 1.0 / 8 + sqrt(30.0) / 144;
const double v2 = 0.5*sqrt((15 - 2*sqrt(30.0)) / 35);
const double v3 = v2 * (1.0/6 - sqrt(30.0)/24);
const double v4 = v2 * (1.0/21 - 5*sqrt(30.0)/168);
const double v5 = v2 - 2*v3;


const std::vector<Vector> GaussLegendre::c_gauss = {
    { (3.0 - std::sqrt(3.0)) / 6.0, (3.0 + std::sqrt(3.0)) / 6.0 },
    { (5.0 - std::sqrt(15.0)) / 10.0, 0.5, (5.0 + std::sqrt(15.0)) / 10.0 },
    { 0.5-w2, 0.5-v2, 0.5+v2, 0.5+w2},
    { 0.0469101, 0.230765, 0.5, 0.769235, 0.95309 }
};

const std::vector<std::vector<Vector>> GaussLegendre::A_gauss = {
    {
        { 0.25, 0.25 - std::sqrt(3.0)/6.0 },
        { 0.25 + std::sqrt(3.0)/6.0, 0.25 }
    },
    {
        { 5.0/36.0, 2.0/9.0 - std::sqrt(15.0)/15.0, 5.0/36.0 - std::sqrt(15.0)/30.0 },
        { 5.0/36.0 + std::sqrt(15.0)/24.0, 2.0/9.0, 5.0/36.0 - std::sqrt(15.0)/24.0 },
        { 5.0/36.0 + std::sqrt(15.0)/30.0, 2.0/9.0 + std::sqrt(15.0)/15.0, 5.0/36.0 }
    },
    {
        { w1, v1-w3+v4, v1-w3-v4, w1-w5 },
        { w1-v3+w4, v1, v1-v5, w1-v3-w4 },
        { w1+v3+w4, v1+v5, v1, w1+v3-w4 },
        { w1+w5, v1+w3+v4, v1+w3-v4, w1 }
    },
    {
        {0.0592317, -0.0195704, 0.0112544, -0.00559379, 0.00158811},
        {0.128151, 0.119657, -0.0245921, 0.0103183, -0.00276899},
        {0.113776, 0.260005, 0.142222, -0.0206903, 0.00468715},
        {0.121232, 0.228996, 0.309037, 0.119657, -0.00968756},
        {0.116875, 0.244908, 0.27319, 0.258885, 0.0592317}
    }
};

const std::vector<Vector> GaussLegendre::b_gauss = {
    { 0.5, 0.5 },
    { 5.0/18.0, 4.0/9.0, 5.0/18.0 },
    { 2*w1, 2*v1, 2*v1, 2*w1},
    { 0.118463, 0.239314, 0.284444, 0.239314, 0.118463 }
};


Vector GaussLegendre::computeStep(const Vector& u_n, double t_n, double dt) {
    int s_local = c_gauss[s-2].size();
    std::vector<Vector> y(s_local);
    std::vector<Vector> k(s_local);

    double epsilon = 1e-15; // 容差
    int max_iterations = 100; // 最大迭代次数

    for (int i = 0; i < s_local; ++i) {
        y[i] = u_n;
        for (int j = 0; j < i; ++j) {
            y[i] += dt * A_gauss[s-2][i][j] * k[j];
        }

        int iterations = 0;
        while (iterations < max_iterations) {
            Vector y_prev = y[i];
            k[i] = F(y_prev, t_n + c_gauss[s-2][i] * dt);
            y[i] = u_n;
            for (int j = 0; j <= i; ++j) {
                y[i] += dt * A_gauss[s-2][i][j] * k[j];
            }
            if ((y[i] - y_prev).empty() || std::abs(Inf_norm(y[i] - y_prev)) < epsilon) {
                break;
            }
            iterations++;
        }

        if (iterations >= max_iterations) {
            std::cout << "Gauss-Legendre: Iteration did not converge within " << max_iterations << " iterations at stage " << i << std::endl;
        }
    }

    Vector y_final = u_n;
    for (int i = 0; i < s_local; ++i) {
        y_final += dt * b_gauss[s-2][i] * k[i];
    }

    return y_final;
}

