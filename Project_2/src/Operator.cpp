#include "Operator.hpp"
#include <stdexcept>
#include <cmath>

int idxx(int n, int i, int j) { return i + (n + 1) * j; }

bool IsEnoughPoint2D(int n, int i, int j) {
    if((i-3)<0 || (i+3)>n || (j+3)>n || (j-3)<0) return false;
    else return true;
}

bool IsEnoughPoint1DX(int n, int i, int j) {
    if(((i-3)>=0 && (i+3)<=n))return true;
    else return false;
}

bool IsEnoughPoint1DY(int n, int i, int j) {
    if(((j+3)<=n && (j-3)>=0)) return true;
    else return false;
}

template <int Dim>
Vector Injection<Dim>::operator()(const Vector& v) const {
    if constexpr (Dim == 1) {
        int n = v.size();
        Vector result((n+1)/2,0.0);
        for(int i=0;i<result.size();i++){
            result[i]=v[2*i];
        }
        return result;
    }
    else if constexpr (Dim == 2){ // 对不规则也不用改
        int n = (int)round(sqrt(v.size())) - 1;
        Vector result((n/2+1)*(n/2+1),0.0);
        for(int i=0;i<=n/2;i++){
            for(int j=0;j<=n/2;j++){
                result[idxx(n/2,i,j)] = v[idxx(n,2*i,2*j)];
            }
        }
        return result;
    } 
    else {
        throw std::runtime_error("Unexpected Dim.");
    }
}

template <int Dim>
Vector LinearInterpolation<Dim>::operator()(const Vector& v) const {
    if constexpr (Dim == 1) {
        int n = v.size();
        Vector result(2*n-1,0.0);
        for(int i=0;i<result.size();i++){
            if(i%2==0) result[i]=v[i/2];
            else result[i] = (v[(i-1)/2] + v[(i+1)/2]) / 2;
        }
        return result;
    } 
    else if constexpr (Dim == 2){
        int n = (int)round(sqrt(v.size())) - 1;
        Vector result((n*2+1)*(n*2+1),0.0);
        for(int i=0;i<=2*n;i++){
            for(int j=0;j<=2*n;j++){
                if(i%2==0 && j%2==0) result[idxx(2*n,i,j)] = v[idxx(n,i/2,j/2)];
                if(i%2==1 && j%2==0) result[idxx(2*n,i,j)] = 0.5*v[idxx(n,(i-1)/2,j/2)] + 0.5*v[idxx(n,(i+1)/2,j/2)];
                if(i%2==0 && j%2==1) result[idxx(2*n,i,j)] = 0.5*v[idxx(n,i/2,(j-1)/2)] + 0.5*v[idxx(n,i/2,(j+1)/2)];
                if(i%2==1 && j%2==1) result[idxx(2*n,i,j)] = 0.25*v[idxx(n,(i-1)/2,(j-1)/2)] + 0.25*v[idxx(n,(i+1)/2,(j+1)/2)] + 0.25*v[idxx(n,(i+1)/2,(j-1)/2)] +0.25*v[idxx(n,(i-1)/2,(j+1)/2)];
            }
        }
        return result;
    }     
    else {
        throw std::runtime_error("Unexpected Dim.");
    }
}

template <int Dim>
Vector FullWeighting<Dim>::operator()(const Vector& v) const {
    if constexpr (Dim == 1) {
        int n = v.size();
        Vector result((n+1)/2,0.0);
        for(int i=1;i<result.size()-1;i++){
            result[i]=0.25*v[2*i-1]+0.5*v[2*i]+0.25*v[2*i+1];
        }
        result[0]=v[0];
        result[(n-1)/2]=v[n-1];
        return result;
    }
    else if constexpr (Dim == 2){
        int n = (int)round(sqrt(v.size())) - 1;
        Vector result((n/2+1)*(n/2+1),0.0);
        for(int i=1;i<n/2;i++) {
            for(int j=1;j<n/2;j++){
                result[idxx(n/2,i,j)] += v[idxx(n,2*i,2*j)]*0.25;
                result[idxx(n/2,i,j)] += (v[idxx(n,2*i-1,2*j-1)] + v[idxx(n,2*i+1,2*j-1)] + v[idxx(n,2*i-1,2*j+1)] + v[idxx(n,2*i+1,2*j+1)]) * 0.125 / 2;
                result[idxx(n/2,i,j)] += (v[idxx(n,2*i-1,2*j)] + v[idxx(n,2*i,2*j-1)] + v[idxx(n,2*i+1,2*j)] + v[idxx(n,2*i,2*j+1)]) * 0.125;
            }
        }
        for(int i=1;i<n/2;i++){
            result[idxx(n/2,i,0)] = 0.5*v[idxx(n,2*i,0)] + 0.25*v[idxx(n,2*i-1,0)] + 0.25*v[idxx(n,2*i+1,0)];
            result[idxx(n/2,i,n/2)] = 0.5*v[idxx(n,2*i,n)] + 0.25*v[idxx(n,2*i-1,n)] + 0.25*v[idxx(n,2*i+1,n)];
            result[idxx(n/2,0,i)] = 0.5*v[idxx(n,0,2*i)] + 0.25*v[idxx(n,0,2*i-1)] + 0.25*v[idxx(n,0,2*i+1)];
            result[idxx(n/2,n/2,i)] = 0.5*v[idxx(n,n,2*i)] + 0.25*v[idxx(n,n,2*i-1)] + 0.25*v[idxx(n,n,2*i+1)];
        }
        result[idxx(n/2,0,0)] = v[idxx(n,0,0)];
        result[idxx(n/2,0,n/2)] = v[idxx(n,0,n)];
        result[idxx(n/2,n/2,0)] = v[idxx(n,n,0)];
        result[idxx(n/2,n/2,n/2)] = v[idxx(n,n,n)];
        return result;
    } 
    else {
        throw std::runtime_error("Unexpected Dim.");
    }
}

template <int Dim>
Vector QuadraticInterpolation<Dim>::operator()(const Vector& v) const {
    if constexpr (Dim == 1) {
        int n = v.size();
        if(n <= 2){
            LinearInterpolation<1> l;
            return l(v);
        }

        Vector result(2*n-1,0.0);
        for(int i=0;i<n;i++){
            result[2*i]=v[i];
        }
        result[1] = 0.375*v[0] + 0.75*v[1] - 0.125*v[2];
        result[2*n-3] = 0.375*v[n-1] + 0.75*v[n-2] - 0.125*v[n-3];
        for(int i=1;i<n-2;i++){
            result[2*i+1] = 9.0/16*(v[i] + v[i+1]) - 1.0/16*(v[i-1] + v[i+2]);
        }
        return result;
    } 
    else if constexpr (Dim == 2) {
        int n = (int)round(sqrt(v.size())) - 1;
        Vector result((n*2+1)*(n*2+1),0.0);
        // 为了点尽量取对称格式，对贴近边界的点选择用线性插值,对内部点取四个方向的二次插值的平均
        for(int i=0;i<=2*n;i++){
            for(int j=0;j<=2*n;j++){
                if(i%2==0 && j%2==0){
                    result[idxx(2*n,i,j)] = v[idxx(n,i/2,j/2)];
                }
                if(i%2==0 && j%2==1){
                    if(IsEnoughPoint1DY(2*n,i,j)) result[idxx(2*n,i,j)] = 9.0/16*(v[idxx(n,i/2,(j-1)/2)] + v[idxx(n,i/2,(j+1)/2)]) -1.0/16*(v[idxx(n,i/2,(j-3)/2)] + v[idxx(n,i/2,(j+3)/2)]);
                    else result[idxx(2*n,i,j)] = 0.5*(v[idxx(n,i/2,(j-1)/2)] + v[idxx(n,i/2,(j+1)/2)]);
                }
                if(i%2==1 && j%2==0){
                    if(IsEnoughPoint1DX(2*n,i,j)) result[idxx(2*n,i,j)] = 9.0/16*(v[idxx(n,(i-1)/2,j/2)] + v[idxx(n,(i+1)/2,j/2)]) -1.0/16*(v[idxx(n,(i-3)/2,j/2)] + v[idxx(n,(i+3)/2,j/2)]);
                    else result[idxx(2*n,i,j)] = 0.5*(v[idxx(n,(i+1)/2,j/2)] + v[idxx(n,(i-1)/2,j/2)]);
                }
                if(i%2==1 && j%2==1){
                    if(IsEnoughPoint2D(2*n,i,j)) result[idxx(2*n,i,j)] = 5.0/16*(v[idxx(n,(i-1)/2,(j-1)/2)] + v[idxx(n,(i-1)/2,(j+1)/2)] + v[idxx(n,(i+1)/2,(j-1)/2)] + v[idxx(n,(i+1)/2,(j+1)/2)]) - 1.0/32*(v[idxx(n,(i-3)/2,(j-1)/2)] + v[idxx(n,(i-3)/2,(j+1)/2)] + v[idxx(n,(i+3)/2,(j-1)/2)] + v[idxx(n,(i+3)/2,(j+1)/2)] + v[idxx(n,(i-1)/2,(j+3)/2)] + v[idxx(n,(i+1)/2,(j+3)/2)] + v[idxx(n,(i+1)/2,(j-3)/2)] + v[idxx(n,(i-1)/2,(j+3)/2)]);
                    else result[idxx(2*n,i,j)] = 0.25*(v[idxx(n,(i-1)/2,(j-1)/2)] + v[idxx(n,(i+1)/2,(j-1)/2)] + v[idxx(n,(i+1)/2,(j+1)/2)] + v[idxx(n,(i-1)/2,(j+1)/2)]);
                }
            }
        }
        return result;
    }
    else {
        throw std::runtime_error("Unexpected Dim.");
    }
}

// 显式实例化模板（仅Dim=1）
template class Injection<1>;
template class LinearInterpolation<1>;
template class FullWeighting<1>;
template class QuadraticInterpolation<1>;
template class Injection<2>;
template class LinearInterpolation<2>;
template class FullWeighting<2>;
template class QuadraticInterpolation<2>;