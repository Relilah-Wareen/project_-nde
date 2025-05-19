#include "../src/IVPFactory.hpp"
#include <fstream>
#include <cmath>
#include <iostream>

Vector rhs_orbit(const Vector& u, double t) {
    const double mu = 0.012277471;
    const double r1_sq = pow(u[0] + mu - 1, 2) + pow(u[1], 2) + pow(u[2], 2);
    const double r2_sq = pow(u[0] + mu, 2) + pow(u[1], 2) + pow(u[2], 2);
    const double r1_cubed = pow(r1_sq, 1.5);
    const double r2_cubed = pow(r2_sq, 1.5);

    return {
        u[3], u[4], u[5],
        2*u[4] + u[0] - (mu*(u[0]+mu-1))/r1_cubed - ((1-mu)*(u[0]+mu))/r2_cubed,
        -2*u[3] + u[1] - (mu*u[1])/r1_cubed - ((1-mu)*u[1])/r2_cubed,
        -(mu*u[2])/r1_cubed - ((1-mu)*u[2])/r2_cubed
    };
}

void run_test(const std::string& method, int p_base, int N_base, double T) {
    auto& factory = IVPFactory::getInstance();
    std::ofstream out("../data/convergence_results.csv", std::ios::app);

    for (int mult : {1, 2, 4}) {
        int N = N_base * mult;
        try {
            auto solver = factory.createSolver(method, p_base, T, rhs_orbit);
            solver->initialize({0.994, 0, 0, 0, -2.0015851063790825224, 0}, N);
            solver->solve();
            
            double error = solver->computePeriodError();
            out << method << "," << p_base << "," << N << "," << error << "\n";
        }
        catch (const std::exception& e) {
            std::cerr << "Error in " << method << " p=" << p_base 
                     << " N=" << N << ": " << e.what() << "\n";
        }
    }
}

int main() {
    RegisterAllSolvers();
    const double T = 17.06521656015796;
    const int N_base = 100000;

    std::ofstream out("../data/convergence_results.csv");
    out << "Method,p,N,Error\n";  // CSV header
    out.close();

    // Adams-Bashforth tests
    for (int p : {1, 2, 3, 4}) {
        run_test("AdamsBashforth", p, N_base, T);
    }

    // Adams-Moulton tests
    for (int p : {2, 3, 4, 5}) {
        run_test("AdamsMoulton", p, N_base, T);
    }

    // BDF tests
    for (int p : {1, 2, 3, 4}) {
        run_test("BDF", p, N_base, T);
    }

    return 0;
}