#ifndef PROBLEMDOMAIN_HPP
#define PROBLEMDOMAIN_HPP

#include <array>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

enum class BCType {
    Dirichlet,
    Neumann,
    Mixed
};

enum class MixedType {
    DDDD, DDDN, DDND, DDNN,
    DNDD, DNDN, DNND, DNNN,
    NDDD, NDDN, NDND, NDNN,
    NNDD, NNDN, NNND, NNNN
};

enum class DomainType {
    Default,  // 默认方形区域
    Circular  // 带有圆形区域的区域
};

class ProblemDomain {
private:
    DomainType domtype;
    int GridNum;
    double circleCenterX, circleCenterY, circleRadius;

    BCType squareBoundaryType;
    BCType circleBoundaryType;
    MixedType squareMixedType;

    std::array<std::function<double(double, double)>, 4> squareBoundaryFunctions;
    std::function<double(double, double)> circleBoundaryFunction;
    std::function<double(double, double)> circleBoundaryFunctionXDerivative; // x方向导数
    std::function<double(double, double)> circleBoundaryFunctionYDerivative; // y方向导数
    std::array<BCType, 4> getEdgeTypes(MixedType type) const; // 四条边的边界类型，顺序为左右下上.
    
    bool isCircleValid() const; // 判断有无四个以上的点在内部.

    // 用于存储 (0,0) 处的值
    double valueAt00 = std::numeric_limits<double>::max();
public:
    ProblemDomain(DomainType type, int gridNum, double centerX = 0.0, double centerY = 0.0, double radius = 0.0);

    void setSquareBoundaryCondition(BCType type, std::array<std::function<double(double, double)>, 4> funcs);
    void setSquareBoundaryCondition(MixedType mixedType, std::array<std::function<double(double, double)>, 4> funcs);
    void setCircleBoundaryCondition(BCType type, std::function<double(double, double)> func);
    void setCircleBoundaryCondition(BCType type, std::function<double(double, double)> funcXDerivative, std::function<double(double, double)> funcYDerivative); // Neumann 边界条件
    
    double getSquareBoundaryValue(int boundaryIndex, double x, double y) const;
    double getCircleBoundaryValue(double x, double y) const;
    
    bool isPointInsideCircle(double x, double y) const;
    void printInfo() const;
    std::array<BCType, 4> getEdgeTypes() const;
    BCType getsquareBoundaryType() const;
    BCType getcircleBoundaryType() const;
    
    DomainType getDomainType() const { return domtype; }
    int getGridNum() const { return GridNum; }

    double getYOnCircle(double x, double minY, double maxY) const;
    double getXOnCircle(double y, double minX, double maxX) const;
    std::array<double, 2> getRayCircleIntersection(double x, double y) const;
    std::array<double, 2> getRaySquareIntersection(int i, int j) const;

    // 设置 (0,0) 处的值
    void SetValueAt00(double value);
    double getValueAt00() const;

    void exportBoundaryPointsToCSV(const std::string& filename) const;
};

#endif // PROBLEMDOMAIN_HPP