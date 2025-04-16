#ifndef BOUNDCON_HPP
#define BOUNDCON_HPP

#include <functional>
#include <variant>
#include <memory>
#include "Function.hpp" // 包含之前实现的Function类
enum class BCType { Dirichlet, Neumann, Mixed};

// 一维边界组合 (左,右)
enum class Mixed1DType { DD, DN, ND, NN, NOTMIX };

// 二维边界组合 (左,右,下,上)
enum class Mixed2DType {
    DDDD, DDDN, DDND, DDNN,
    DNDD, DNDN, DNND, DNNN,
    NDDD, NDDN, NDND, NDNN,
    NNDD, NNDN, NNND, NNNN,
    NOTMIX
};

class BoundCon {
    private:
        int dimension;
        BCType bc_type;
        Mixed1DType mix1d;
        Mixed2DType mix2d;
    
        // 1D存储（固定值）
        double left_val_1d;
        double right_val_1d;
    
        // 2D直接存储Function对象
        Function left_func_2d;
        Function right_func_2d;
        Function bottom_func_2d;
        Function top_func_2d;
        BCType left_tp;
        BCType right_tp;
        BCType bottom_tp;
        BCType top_tp;
        bool IrregularFlag;
    public:
        // 1D构造函数
        BoundCon(int dim, BCType type, double left, double right, 
                 Mixed1DType m1d = Mixed1DType::NOTMIX);
        
        // 2D Dirichlet构造函数
        BoundCon(int dim, BCType type, const Function& uniform_func);
        
        // 2D Neumann构造函数
        BoundCon(int dim, BCType type,
            const Function& left, const Function& right,
            const Function& bottom, const Function& top);

        // 2D Mixed构造函数
        BoundCon(int dim, BCType type, Mixed2DType m2d,
                const Function& left, const Function& right,
                const Function& bottom, const Function& top);
    
        // 获取方法
        double getLeft1D() const;
        double getRight1D() const;
        double getBound2D(double x, double y) const;
        int getDim() const { return dimension; }
        BCType getBCType() const { return bc_type; }
        BCType getLeft2DType() const { return left_tp; }
        BCType getRight2DType() const { return right_tp; }
        BCType getBottom2DType() const { return bottom_tp; }
        BCType getTop2DType() const { return top_tp; }
        void set_Irregular();
        bool IsIrregular() const { return IrregularFlag; }

        Mixed1DType getMixed1DType() const { 
            if (dimension != 1) throw std::runtime_error("Not a 1D boundary");
            return mix1d; 
        }
        Mixed2DType getMixed2DType() const { 
            if (dimension != 2) throw std::runtime_error("Not a 2D boundary");
            return mix2d; 
        }
    };
    
#endif // BOUNDCON_HPP