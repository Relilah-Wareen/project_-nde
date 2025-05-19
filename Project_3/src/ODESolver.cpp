#include "ODESolver.hpp"
#include <algorithm>
#include <cmath>

ODESolver::ODESolver(int p, double T, std::function<Vector(const Vector&, double)> f) 
    : p(p), T(T), F(f) {}

const std::vector<Vector>& ODESolver::get_results() const {
    return result;
}

void ODESolver::printforU1U2(const std::string& filename) const {
    if (Dim < 2) throw std::runtime_error("Solution dimension must be >= 2");
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open file: " + filename);
    out << "u1,u2\n";
    for (int i = 0; i < result.size(); ++i) {
        out << result[i][0] << "," << result[i][1] << "\n";
    }
}

double ODESolver::computePeriodError(){
    if(!result.size()) return 0;
    else {
        int resz = result.size();
        Vector r = result[0] - result[resz-1];
        return L2_norm(r) / L2_norm(result[0]);
    }
}