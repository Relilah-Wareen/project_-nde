#ifndef EMBEDDED_RK_METHOD_HPP
#define EMBEDDED_RK_METHOD_HPP

#include "ODESolver.hpp"
#include <vector>
#include <cmath>

class EmbeddedRKMethod : public ODESolver {
protected:
    double Ea;                    
    double Er;                    
    double rho_max;              
    double rho_min;               
    double rho;                   
    double current_dt;           

    virtual int getP1() const = 0;          
    virtual int getP2() const = 0;          
    virtual const std::vector<Vector>& getA() const = 0;  
    virtual const Vector& getC() const = 0;               
    virtual const Vector& getB1() const = 0;            
    virtual const Vector& getB2() const = 0;             

public:
    EmbeddedRKMethod(double T, std::function<Vector(const Vector&, double)> f)
        : ODESolver(0, T, f),  // 基类阶数设为0，子类实际阶数通过getP1()决定
          Ea(1e-5), Er(1e-7), rho_max(5.0), rho_min(0.2), rho(0.9), current_dt(0.0) {}

    void initialize(const Vector& init_value,int n=0) override {
        Dim = init_value.size();
        result.resize(0, Vector(Dim, 0.0));
        result.push_back(init_value);
        current_dt = 0.01;
    }

    void solve() override;
};

class Fehlberg45 : public EmbeddedRKMethod {
    public:
        Fehlberg45(double T, std::function<Vector(const Vector&, double)> f)
            : EmbeddedRKMethod(T, f) {}
    
    protected:
        int getP1() const override { return 4; }  // 高阶方法阶数
        int getP2() const override { return 5; }  // 低阶方法阶数
    
        // Butcher表系数（Fehlberg 4(5)）
        const std::vector<Vector>& getA() const override {
            static const std::vector<Vector> A = {
                {},  
                {1.0/4.0},
                {3.0/32.0, 9.0/32.0},
                {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0},
                {439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0},
                {-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0}
            };
            return A;
        }
    
        const Vector& getC() const override {
            static const Vector c = {0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 0.5};
            return c;
        }
    
        const Vector& getB1() const override {  // 4阶方法权重
            static const Vector b1 = {
                25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -0.2, 0.0
            };
            return b1;
        }
    
        const Vector& getB2() const override {  // 5阶方法权重
            static const Vector b2 = {
                16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0
            };
            return b2;
        }
    };
    
    
    class DormandPrince54 : public EmbeddedRKMethod {
    public:
        DormandPrince54(double T, std::function<Vector(const Vector&, double)> f)
            : EmbeddedRKMethod(T, f) {}
    
    protected:
        int getP1() const override { return 5; }  // 高阶方法阶数
        int getP2() const override { return 4; }  // 低阶方法阶数
    
        // Butcher表系数（Dormand-Prince 5(4)）
        const std::vector<Vector>& getA() const override {
            static const std::vector<Vector> A = {
                {}, 
                {1.0/5.0},
                {3.0/40.0, 9.0/40.0},
                {44.0/45.0, -56.0/15.0, 32.0/9.0},
                {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0},
                {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0},
                {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0}
            };
            return A;
        }
    
        const Vector& getC() const override {
            static const Vector c = {0.0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0};
            return c;
        }
    
        const Vector& getB1() const override {  // 5阶方法权重
            static const Vector b1 = {
                35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0
            };
            return b1;
        }
    
        const Vector& getB2() const override {  // 4阶方法权重
            static const Vector b2 = {
                5179.0/57600.0, 0.0, 7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0
            };
            return b2;
        }
    };


#endif // EMBEDDED_RK_METHOD_HPP