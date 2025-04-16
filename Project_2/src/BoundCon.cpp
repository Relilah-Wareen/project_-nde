#include "BoundCon.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
// 1D构造函数
BoundCon::BoundCon(int dim, BCType type, double left, double right, Mixed1DType m1d)
    : dimension(dim), bc_type(type), mix1d(m1d), mix2d(Mixed2DType::NOTMIX),
      left_val_1d(left), right_val_1d(right) 
{
    if (dim != 1) throw std::invalid_argument("Dimension must be 1");
    if (type == BCType::Mixed && m1d == Mixed1DType::NOTMIX) {
        throw std::invalid_argument("Must specify Mixed1DType for mixed BC");
    }
    if(mix1d==Mixed1DType::DD) bc_type = BCType::Dirichlet;
    else if(mix1d==Mixed1DType::NN) bc_type = BCType::Neumann;
    else {}
    IrregularFlag = false;
}

// 2D Dirichlet构造函数
BoundCon::BoundCon(int dim, BCType type, const Function& uniform_func)
    : dimension(dim), bc_type(type), mix1d(Mixed1DType::NOTMIX), mix2d(Mixed2DType::NOTMIX),
      left_func_2d(uniform_func),
      right_func_2d(uniform_func),
      bottom_func_2d(uniform_func),
      top_func_2d(uniform_func),
      left_tp(type),
      right_tp(type),
      bottom_tp(type),
      top_tp(type) 
{
    if (dim != 2) throw std::invalid_argument("Dimension must be 2");
    if (type != BCType::Dirichlet) throw std::invalid_argument("Use Neumann & Mixed constructor");
    IrregularFlag = false;
}

BoundCon::BoundCon(int dim, BCType type,
    const Function& left, const Function& right,
    const Function& bottom, const Function& top)
    : dimension(dim), bc_type(type), mix1d(Mixed1DType::NOTMIX), mix2d(Mixed2DType::NOTMIX),
    left_func_2d(left),
    right_func_2d(right),
    bottom_func_2d(bottom),
    top_func_2d(top),
    left_tp(type),
    right_tp(type),
    bottom_tp(type),
    top_tp(type) 
{
    if (dim != 2) throw std::invalid_argument("Dimension must be 2");
    if (type != BCType::Neumann) throw std::invalid_argument("Use Dirichlet & Mixed constructor");
    IrregularFlag = false;
}
// 2D Mixed构造函数
BoundCon::BoundCon(int dim, BCType type, Mixed2DType m2d,
                  const Function& left, const Function& right,
                  const Function& bottom, const Function& top)
    : dimension(dim), bc_type(type), mix1d(Mixed1DType::NOTMIX), mix2d(m2d),
      left_func_2d(left),
      right_func_2d(right),
      bottom_func_2d(bottom),
      top_func_2d(top) 
{
    if (dim != 2) throw std::invalid_argument("Dimension must be 2");
    if (type != BCType::Mixed) throw std::invalid_argument("This is Mixed constructor");
    switch (m2d) {
        case Mixed2DType::DDDD: {left_tp=BCType::Dirichlet, right_tp=BCType::Dirichlet, bottom_tp=BCType::Dirichlet, top_tp=BCType::Dirichlet;bc_type=BCType::Dirichlet;break;}
        case Mixed2DType::DDDN: {left_tp=BCType::Dirichlet, right_tp=BCType::Dirichlet, bottom_tp=BCType::Dirichlet, top_tp=BCType::Neumann;break;}
        case Mixed2DType::DDND: {left_tp=BCType::Dirichlet, right_tp=BCType::Dirichlet, bottom_tp=BCType::Neumann, top_tp=BCType::Dirichlet;break;}
        case Mixed2DType::DDNN: {left_tp=BCType::Dirichlet, right_tp=BCType::Dirichlet, bottom_tp=BCType::Neumann, top_tp=BCType::Neumann;break;}
        case Mixed2DType::DNDD: {left_tp=BCType::Dirichlet, right_tp=BCType::Neumann, bottom_tp=BCType::Dirichlet, top_tp=BCType::Dirichlet;break;}
        case Mixed2DType::DNDN: {left_tp=BCType::Dirichlet, right_tp=BCType::Neumann, bottom_tp=BCType::Dirichlet, top_tp=BCType::Neumann;break;}
        case Mixed2DType::DNND: {left_tp=BCType::Dirichlet, right_tp=BCType::Neumann, bottom_tp=BCType::Neumann, top_tp=BCType::Dirichlet;break;}
        case Mixed2DType::DNNN: {left_tp=BCType::Dirichlet, right_tp=BCType::Neumann, bottom_tp=BCType::Neumann, top_tp=BCType::Neumann;break;}
        case Mixed2DType::NDDD: {left_tp=BCType::Neumann, right_tp=BCType::Dirichlet, bottom_tp=BCType::Dirichlet, top_tp=BCType::Dirichlet;break;}
        case Mixed2DType::NDDN: {left_tp=BCType::Neumann, right_tp=BCType::Dirichlet, bottom_tp=BCType::Dirichlet, top_tp=BCType::Neumann;break;}
        case Mixed2DType::NDND: {left_tp=BCType::Neumann, right_tp=BCType::Dirichlet, bottom_tp=BCType::Neumann, top_tp=BCType::Dirichlet;break;}
        case Mixed2DType::NDNN: {left_tp=BCType::Neumann, right_tp=BCType::Dirichlet, bottom_tp=BCType::Neumann, top_tp=BCType::Neumann;break;}
        case Mixed2DType::NNDD: {left_tp=BCType::Neumann, right_tp=BCType::Neumann, bottom_tp=BCType::Dirichlet, top_tp=BCType::Dirichlet;break;}
        case Mixed2DType::NNDN: {left_tp=BCType::Neumann, right_tp=BCType::Neumann, bottom_tp=BCType::Dirichlet, top_tp=BCType::Neumann;break;}
        case Mixed2DType::NNND: {left_tp=BCType::Neumann, right_tp=BCType::Neumann, bottom_tp=BCType::Neumann, top_tp=BCType::Dirichlet;break;}
        case Mixed2DType::NNNN: {left_tp=BCType::Neumann, right_tp=BCType::Neumann, bottom_tp=BCType::Neumann, top_tp=BCType::Neumann;bc_type=BCType::Neumann;break;}
        default: throw std::invalid_argument("Unknown MixedType");
    }
    IrregularFlag = false;
}

// 获取边界值方法
double BoundCon::getLeft1D() const {
    if (dimension != 1) throw std::runtime_error("Not a 1D boundary");
    return left_val_1d;
}

double BoundCon::getRight1D() const {
    if (dimension != 1) throw std::runtime_error("Not a 1D boundary");
    return right_val_1d;
}

double BoundCon::getBound2D(double x, double y) const {
    if (dimension != 2) throw std::runtime_error("Not a 2D boundary");
    if(IrregularFlag){
        const double eps = 1e-10;
        bool on_left   = (x <= eps);
        bool on_right  = (x >= 1.0 - eps);
        bool on_bottom = (y >= 1.0/16*sin(M_PI*x)-eps && y<=1.0/16*sin(M_PI*x)+eps );
        bool on_top    = (y >= 1.0 - eps);

        if (on_bottom) return bottom_func_2d(x, y);
        if (on_top)    return top_func_2d(x, 1.0);
        if (on_left)   return left_func_2d(0.0, y);
        if (on_right)  return right_func_2d(1.0, y);
    }   
    else {
        const double eps = 1e-10;
        bool on_left   = (x <= eps);
        bool on_right  = (x >= 1.0 - eps);
        bool on_bottom = (y <= eps);
        bool on_top    = (y >= 1.0 - eps);

        if (on_bottom) return bottom_func_2d(x, 0.0);
        if (on_top)    return top_func_2d(x, 1.0);    
        if (on_left)   return left_func_2d(0.0, y);
        if (on_right)  return right_func_2d(1.0, y);

        throw std::runtime_error("Point is not on boundary");
    } 
    return 0;
}

void BoundCon::set_Irregular() {
    if(dimension != 2) throw std::runtime_error("Unexpected Dim. Irregular problem only for 2D in this project.");
    else IrregularFlag = true; 
}
