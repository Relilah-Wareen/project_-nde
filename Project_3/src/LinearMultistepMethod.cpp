#include "LinearMultistepMethod.hpp"

LinearMultistepMethod::LinearMultistepMethod(int p, double T, std::function<Vector(const Vector&, double)> f) 
    : ODESolver(p, T, f) {} 

// RK4单步实现
Vector LinearMultistepMethod::RK_4onestep(const Vector& u_n, double t) {
    Vector y1 = F(u_n, t);
    Vector y2 = F(u_n + 0.5 * k * y1, t + 0.5 * k);
    Vector y3 = F(u_n + 0.5 * k * y2, t + 0.5 * k);
    Vector y4 = F(u_n + k * y3, t + k);
    return u_n + (k / 6.0) * (y1 + 2*y2 + 2*y3 + y4);
}

// 初始化方法
void LinearMultistepMethod::initialize(const Vector& init_value, int n) {
    Dim = init_value.size();
    Num = n;
    k = T / n;
    result.resize(n + 1, Vector(Dim, 0.0));
    result[0] = init_value;
    for (int i = 1; i < s; ++i) {  // 生成s个初值
        result[i] = RK_4onestep(result[i-1], (i-1)*k);
    }
}

AdamsMethod::AdamsMethod(const int& p, const double& T, const std::function<Vector(const Vector&, double)>& f)
    : LinearMultistepMethod(p, T, f) {}

AdamsBashforth::AdamsBashforth(const int& p, const double& T, const std::function<Vector(const Vector&, double)>& f)
    : AdamsMethod(p, T, f) {
        s = p; 
        if(p > 4 || p < 1) throw std::invalid_argument("p = 1,2,3,4.");
    }
    
const std::unordered_map<int, std::vector<double>> AdamsBashforth::beta_coeffs = {
    {1, {1.0}},
    {2, {-1.0/2, 3.0/2}},
    {3, {5.0/12, -16.0/12, 23.0/12}},
    {4, {-9.0/24, 37.0/24, -59.0/24, 55.0/24}}
};

double AdamsBashforth::beta(int i, int j) const {
    try {
        return beta_coeffs.at(i).at(j);
    } 
    catch (const std::out_of_range& e) {
        throw std::out_of_range("Invalid beta coefficient index: i=" + 
                              std::to_string(i) + " j=" + std::to_string(j));
    }
}

void AdamsBashforth::solve() {
    for(int i = s; i <= Num; i++){
        Vector rhs(Dim,0.0);
        int idx = i-s;
        for(int j=0; j<s; j++){
            rhs += beta(s,j) * F(result[idx+j], (idx+j) * k);
        }
        result[i] = result[i-1] + 1.0*k*rhs;
    }
} 

AdamsMoulton::AdamsMoulton(const int& p, const double& T, const std::function<Vector(const Vector&, double)>& f)
    : AdamsMethod(p, T, f) {
        s = p - 1;
        if(p > 5 || p < 2) throw std::invalid_argument("p = 2,3,4,5.");
    }

const std::unordered_map<int, std::vector<double>> AdamsMoulton::beta_coeffs = {
    {1, {1.0/2, 1.0/2}},
    {2, {-1.0/12, 8.0/12, 5.0/12}},
    {3, {1.0/24, -5.0/24, 19.0/24, 9.0/24}},
    {4, {-19.0/720, 106.0/720, -264.0/720, 646.0/720, 251.0/720}}
};

double AdamsMoulton::beta(int i, int j) const {
    try {
        return beta_coeffs.at(i).at(j);
    } 
    catch (const std::out_of_range& e) {
        throw std::out_of_range("Invalid beta coefficient index: i=" + 
                              std::to_string(i) + " j=" + std::to_string(j));
    }
}

/*
对本次项目情况中的隐式方法，我们目前只实现了不动点迭代的方法，对于可能不收敛的情况，会输出可能异常的提示。
*/
void AdamsMoulton::solve() {
    double epsilon = 1e-15;
    for(int i = s; i <= Num; i++){
        Vector u = result[i-1];
        int idx = i-s;
        int cnt = 0;
        while (cnt<50){
            Vector curr(Dim,0.0);
            Vector rhs(Dim,0.0);
            for(int j=0; j<s; j++){
                rhs += beta(s,j) * F(result[idx+j], (idx+j) * k);
            }
            rhs += beta(s,s) * F(u, i * k);
            curr = result[i-1] + k * rhs;
            if(std::abs(Inf_norm(curr - u)) < epsilon) {
                u=curr;cnt++;
                break;
            }
            cnt++;
            u=curr;
        }
        if(cnt>=50)std::cout << "Implicit method iterative solution may not converge." << std::endl;
        result[i] = u;
    }
}

BDFMethod::BDFMethod(const int& p, const double& T, const std::function<Vector(const Vector&, double)>& f)
    : LinearMultistepMethod(p, T, f) {
    s = p; 
    if(p > 4 || p < 1) throw std::invalid_argument("p = 1,2,3,4.");            
}

const std::unordered_map<int, Vector> BDFMethod::alpha_coeffs = {
    {1, {-1.0, 1.0}},
    {2, {1.0/3, -4.0/3, 1.0}},
    {3, {-2.0/11, 9.0/11, -18.0/11, 1.0}},
    {4, {3.0/25, -16.0/25, 36.0/25, -48.0/25, 1.0}}
};

const std::unordered_map<int, double> BDFMethod::beta_s = {
    {1, 1.0},
    {2, 2.0/3},
    {3, 6.0/11},
    {4, 12.0/25}
};

double BDFMethod::alpha(int i, int j) {
    try {
        return alpha_coeffs.at(i).at(j);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range("Invalid alpha coefficient index: i=" + std::to_string(i) + " j=" + std::to_string(j));
    }
}

double BDFMethod::beta(int s) {
    try {
        return beta_s.at(s);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range("Invalid beta index: s=" + std::to_string(s));
    }
}

void BDFMethod::solve() {
    double epsilon = 1e-15;
    for(int i = s; i <= Num; i++){
        Vector u = result[i-1];
        int idx = i-s;
        int cnt = 0;
        while (cnt<50){
            Vector curr(Dim,0.0);
            // Vector rhs(Dim,0.0);
            for(int j=0; j<s; j++){
                curr += -1.0 * alpha(s,j) * result[idx + j];
            }
            curr += 1.0 * k * beta(s) * F(u, i * k);
            // curr = rhs;
            if(std::abs(Inf_norm(curr - u)) < epsilon) {
                u=curr;cnt++;
                break;
            }
            cnt++;
            u=curr;
        }
        if(cnt>=50)std::cout << "Implicit method iterative solution may not converge." << std::endl;
        result[i] = u;
    }
}