#include "EmbeddedRKMethod.hpp"
#include <algorithm>

void EmbeddedRKMethod::solve() {
    if (result.empty()) {
        throw std::runtime_error("Call initialize() before solve().");
    }

    double t_current = 0.0;
    size_t step = 0;
    Vector u_prev = result[0];

    while (t_current < T) {
        // (a) 用两个方法计算u1和u2
        const auto& A = getA();
        const auto& c = getC();
        const auto& b1 = getB1();
        const auto& b2 = getB2();

        std::vector<Vector> k(A.size(), Vector(Dim, 0.0));
        for (size_t i = 0; i < A.size(); ++i) {
            Vector sum = u_prev;
            for (size_t j = 0; j < i; ++j) {
                sum += current_dt * A[i][j] * k[j];
            }
            k[i] = F(sum, t_current + c[i] * current_dt);
        }

        Vector u1 = u_prev;
        Vector u2 = u_prev;
        for (size_t i = 0; i < A.size(); ++i) {
            u1 += current_dt * b1[i] * k[i];
            u2 += current_dt * b2[i] * k[i];
        }

        // (b) 计算误差e
        Vector epsilon(Dim);
        for (size_t i = 0; i < Dim; ++i) {
            epsilon[i] = Ea + std::abs(u_prev[i]) * Er;
        }
        double e_sq = 0.0;
        for (size_t i = 0; i < Dim; ++i) {
            double diff = (u1[i] - u2[i]) / epsilon[i];
            e_sq += diff * diff;
        }
        double e = std::sqrt(e_sq / Dim);

        // (c) 调整步长dt
        int q = std::min(getP1(), getP2());
        double factor = rho * std::pow(1.0 / e, 1.0 / (q + 1));
        factor = std::max(rho_min, std::min(rho_max, factor));
        double new_dt = current_dt * factor;
        new_dt = std::min(new_dt, T - t_current);  // 确保不超出T

        // (d) 判断是否接受当前步
        if (e <= 1.0) {
            result.push_back(u1);
            u_prev = u1;
            t_current += current_dt;
            current_dt = new_dt;
            step++;
        } else {
            current_dt = new_dt;
        }
    }
}


