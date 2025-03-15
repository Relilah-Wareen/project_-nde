#include "ProblemDomain.hpp"

std::array<BCType, 4> ProblemDomain::getEdgeTypes(MixedType type) const {
    switch (type) {
        case MixedType::DDDD: return {BCType::Dirichlet, BCType::Dirichlet, BCType::Dirichlet, BCType::Dirichlet};
        case MixedType::DDDN: return {BCType::Dirichlet, BCType::Dirichlet, BCType::Dirichlet, BCType::Neumann};
        case MixedType::DDND: return {BCType::Dirichlet, BCType::Dirichlet, BCType::Neumann, BCType::Dirichlet};
        case MixedType::DDNN: return {BCType::Dirichlet, BCType::Dirichlet, BCType::Neumann, BCType::Neumann};
        case MixedType::DNDD: return {BCType::Dirichlet, BCType::Neumann, BCType::Dirichlet, BCType::Dirichlet};
        case MixedType::DNDN: return {BCType::Dirichlet, BCType::Neumann, BCType::Dirichlet, BCType::Neumann};
        case MixedType::DNND: return {BCType::Dirichlet, BCType::Neumann, BCType::Neumann, BCType::Dirichlet};
        case MixedType::DNNN: return {BCType::Dirichlet, BCType::Neumann, BCType::Neumann, BCType::Neumann};
        case MixedType::NDDD: return {BCType::Neumann, BCType::Dirichlet, BCType::Dirichlet, BCType::Dirichlet};
        case MixedType::NDDN: return {BCType::Neumann, BCType::Dirichlet, BCType::Dirichlet, BCType::Neumann};
        case MixedType::NDND: return {BCType::Neumann, BCType::Dirichlet, BCType::Neumann, BCType::Dirichlet};
        case MixedType::NDNN: return {BCType::Neumann, BCType::Dirichlet, BCType::Neumann, BCType::Neumann};
        case MixedType::NNDD: return {BCType::Neumann, BCType::Neumann, BCType::Dirichlet, BCType::Dirichlet};
        case MixedType::NNDN: return {BCType::Neumann, BCType::Neumann, BCType::Dirichlet, BCType::Neumann};
        case MixedType::NNND: return {BCType::Neumann, BCType::Neumann, BCType::Neumann, BCType::Dirichlet};
        case MixedType::NNNN: return {BCType::Neumann, BCType::Neumann, BCType::Neumann, BCType::Neumann};
        default: throw std::invalid_argument("Unknown MixedType");
    }
}

std::array<BCType, 4> ProblemDomain::getEdgeTypes() const {
    switch (squareMixedType) {
        case MixedType::DDDD: return {BCType::Dirichlet, BCType::Dirichlet, BCType::Dirichlet, BCType::Dirichlet};
        case MixedType::DDDN: return {BCType::Dirichlet, BCType::Dirichlet, BCType::Dirichlet, BCType::Neumann};
        case MixedType::DDND: return {BCType::Dirichlet, BCType::Dirichlet, BCType::Neumann, BCType::Dirichlet};
        case MixedType::DDNN: return {BCType::Dirichlet, BCType::Dirichlet, BCType::Neumann, BCType::Neumann};
        case MixedType::DNDD: return {BCType::Dirichlet, BCType::Neumann, BCType::Dirichlet, BCType::Dirichlet};
        case MixedType::DNDN: return {BCType::Dirichlet, BCType::Neumann, BCType::Dirichlet, BCType::Neumann};
        case MixedType::DNND: return {BCType::Dirichlet, BCType::Neumann, BCType::Neumann, BCType::Dirichlet};
        case MixedType::DNNN: return {BCType::Dirichlet, BCType::Neumann, BCType::Neumann, BCType::Neumann};
        case MixedType::NDDD: return {BCType::Neumann, BCType::Dirichlet, BCType::Dirichlet, BCType::Dirichlet};
        case MixedType::NDDN: return {BCType::Neumann, BCType::Dirichlet, BCType::Dirichlet, BCType::Neumann};
        case MixedType::NDND: return {BCType::Neumann, BCType::Dirichlet, BCType::Neumann, BCType::Dirichlet};
        case MixedType::NDNN: return {BCType::Neumann, BCType::Dirichlet, BCType::Neumann, BCType::Neumann};
        case MixedType::NNDD: return {BCType::Neumann, BCType::Neumann, BCType::Dirichlet, BCType::Dirichlet};
        case MixedType::NNDN: return {BCType::Neumann, BCType::Neumann, BCType::Dirichlet, BCType::Neumann};
        case MixedType::NNND: return {BCType::Neumann, BCType::Neumann, BCType::Neumann, BCType::Dirichlet};
        case MixedType::NNNN: return {BCType::Neumann, BCType::Neumann, BCType::Neumann, BCType::Neumann};
        default: throw std::invalid_argument("Unknown MixedType");
    }
}

ProblemDomain::ProblemDomain(DomainType type, int gridNum, double centerX, double centerY, double radius)
    : domtype(type), GridNum(gridNum), circleCenterX(centerX), circleCenterY(centerY), circleRadius(radius) {
    if (domtype == DomainType::Circular && !isCircleValid()) {
        throw std::invalid_argument("The circular region does not meet the requirements. Please check.");
    }
}

bool ProblemDomain::isCircleValid() const {
    if (domtype == DomainType::Default) return true;

    // 计算 h
    double h = 1.0 / GridNum;

    // 需要数值上在 [h, 1-h] * [h, 1-h] 内部
    double leftBoundary = circleCenterX - circleRadius;
    double rightBoundary = circleCenterX + circleRadius;
    double bottomBoundary = circleCenterY - circleRadius;
    double topBoundary = circleCenterY + circleRadius;

    if (leftBoundary <= h || rightBoundary >= 1 - h ||
        bottomBoundary <= h || topBoundary >= 1 - h) {
        return false;
    }

    // 检查圆形区域内是否至少包含 4 个离散点
    double dx = h;
    double dy = h;
    int pointsInside = 0;

    for (int i = 0; i < GridNum; ++i) {
        for (int j = 0; j < GridNum; ++j) {
            double x = i * dx;
            double y = j * dy;
            if (isPointInsideCircle(x, y)) {
                pointsInside++;
                if (pointsInside >= 4) {
                    return true;
                }
            }
        }
    }

    return false;
}

void ProblemDomain::setSquareBoundaryCondition(BCType type, std::array<std::function<double(double, double)>, 4> funcs) {
    if (type == BCType::Mixed) {
        throw std::invalid_argument("Use setSquareBoundaryCondition with MixedType for mixed boundary conditions.");
    }
    squareBoundaryType = type;
    squareBoundaryFunctions = funcs;
}

void ProblemDomain::setSquareBoundaryCondition(MixedType mixedType, std::array<std::function<double(double, double)>, 4> funcs) {
    if (mixedType == MixedType::DDDD) {
        squareBoundaryType = BCType::Dirichlet;
    } else if (mixedType == MixedType::NNNN) {
        squareBoundaryType = BCType::Neumann;
    } else {
        squareBoundaryType = BCType::Mixed;
        squareMixedType = mixedType;
    }

    squareBoundaryFunctions = funcs;
}

void ProblemDomain::setCircleBoundaryCondition(BCType type, std::function<double(double, double)> func) {
    if (domtype != DomainType::Circular) {
        throw std::invalid_argument("Circle boundary condition can only be set for Circular domain type.");
    }

    if (type != BCType::Dirichlet) {
        throw std::invalid_argument("This overload is for Dirichlet boundary condition only.");
    }

    circleBoundaryType = type;
    circleBoundaryFunction = func;
    circleBoundaryFunctionXDerivative = nullptr;
    circleBoundaryFunctionYDerivative = nullptr; 
}

void ProblemDomain::setCircleBoundaryCondition(BCType type, std::function<double(double, double)> funcXDerivative, std::function<double(double, double)> funcYDerivative) {
    if (domtype != DomainType::Circular) {
        throw std::invalid_argument("Circle boundary condition can only be set for Circular domain type.");
    }

    if (type != BCType::Neumann) {
        throw std::invalid_argument("This overload is for Neumann boundary condition only.");
    }

    circleBoundaryType = type;
    circleBoundaryFunction = nullptr; // 清空Dirichlet函数
    circleBoundaryFunctionXDerivative = funcXDerivative;
    circleBoundaryFunctionYDerivative = funcYDerivative;
}

double ProblemDomain::getSquareBoundaryValue(int boundaryIndex, double x, double y) const {
    // 检查边界索引是否有效
    if (boundaryIndex < 0 || boundaryIndex >= 4) {
        throw std::out_of_range("Invalid boundary index.");
    }

    // 检查 (x, y) 是否位于方形区域的边界上
    bool isOnBoundary = false;
    switch (boundaryIndex) {
        case 0: // 左边界 (x = 0)
            isOnBoundary = (x == 0.0 && y >= 0.0 && y <= 1.0);
            break;
        case 1: // 右边界 (x = 1)
            isOnBoundary = (x == 1.0 && y >= 0.0 && y <= 1.0);
            break;
        case 2: // 下边界 (y = 0)
            isOnBoundary = (y == 0.0 && x >= 0.0 && x <= 1.0);
            break;
        case 3: // 上边界 (y = 1)
            isOnBoundary = (y == 1.0 && x >= 0.0 && x <= 1.0);
            break;
    }

    if (!isOnBoundary) {
        throw std::invalid_argument("The point (x, y) is not on the specified boundary.");
    }

    // 返回边界值
    return squareBoundaryFunctions[boundaryIndex](x, y);
}

double ProblemDomain::getCircleBoundaryValue(double x, double y) const {
    if (domtype != DomainType::Circular) {
        throw std::invalid_argument("Circle boundary condition can only be used for Circular domain type.");
    }

    double dx = x - circleCenterX;
    double dy = y - circleCenterY;
    double distanceSquared = dx * dx + dy * dy;
    double radiusSquared = circleRadius * circleRadius;

    const double epsilon = 1e-6;
    if (std::abs(distanceSquared - radiusSquared) > epsilon) {
        throw std::invalid_argument("The point (x, y) is not on the circular boundary.");
    }

    if (circleBoundaryType == BCType::Dirichlet) {
        return circleBoundaryFunction(x, y);
    } else if (circleBoundaryType == BCType::Neumann) {
        // 计算法向量（单位向量）
        double distance = std::sqrt(distanceSquared);
        double nx = dx / distance; // 法向量的 x 分量
        double ny = dy / distance; // 法向量的 y 分量

        // 获取 x 和 y 方向的偏导数
        double du_dx = circleBoundaryFunctionXDerivative(x, y);
        double du_dy = circleBoundaryFunctionYDerivative(x, y);

        // 计算法向导数：梯度向量与法向量的点积
        double normal_derivative = du_dx * nx + du_dy * ny;
        return normal_derivative;
    } else {
        throw std::invalid_argument("Unexpected circular boundary type.");
    }
}

bool ProblemDomain::isPointInsideCircle(double x, double y) const {
    if (domtype == DomainType::Default) return false;

    double dx = x - circleCenterX;
    double dy = y - circleCenterY;
    return (dx * dx + dy * dy) < (circleRadius * circleRadius);
}

void ProblemDomain::printInfo() const {
    std::cout << "Domain Type: " << (domtype == DomainType::Default ? "Default" : "Circular") << std::endl;
    std::cout << "Grid Number: " << GridNum << std::endl;
    if (domtype == DomainType::Circular) {
        std::cout << "Circle Center: (" << circleCenterX << ", " << circleCenterY << ")" << std::endl;
        std::cout << "Circle Radius: " << circleRadius << std::endl;
    }

    std::cout << "Square Boundary Condition: ";
    if (squareBoundaryType == BCType::Dirichlet) {
        std::cout << "Pure Dirichlet (DDDD)" << std::endl;
    } else if (squareBoundaryType == BCType::Neumann) {
        std::cout << "Pure Neumann (NNNN)" << std::endl;
    } else {
        std::cout << "Mixed Type: ";
        auto edges = getEdgeTypes(squareMixedType);
        for (auto t : edges) {
            std::cout << (t == BCType::Dirichlet ? "D" : "N");
        }
        std::cout << std::endl;
    }

    std::cout << "Circle Boundary Condition: ";
    if (circleBoundaryType == BCType::Dirichlet) {
        std::cout << "Dirichlet" << std::endl;
    } else if (circleBoundaryType == BCType::Neumann) {
        std::cout << "Neumann" << std::endl;
    } else {
        std::cout << "Mixed (Not Supported)" << std::endl;
    }
    
}

BCType ProblemDomain::getsquareBoundaryType() const {
    if (squareBoundaryFunctions[0] == nullptr || 
        squareBoundaryFunctions[1] == nullptr || 
        squareBoundaryFunctions[2] == nullptr || 
        squareBoundaryFunctions[3] == nullptr) {
        throw std::runtime_error("Square boundary condition has not been set.");
    }

    return squareBoundaryType;
}

BCType ProblemDomain::getcircleBoundaryType() const {
    if (circleBoundaryFunction == nullptr && circleBoundaryFunctionXDerivative == nullptr && circleBoundaryFunctionYDerivative == nullptr) {
        throw std::runtime_error("Circle boundary condition has not been set.");
    }

    return circleBoundaryType;
}

    // 根据 x 和 y 的范围，获取圆上的点
double ProblemDomain::getYOnCircle(double x, double minY, double maxY) const {
    if (domtype != DomainType::Circular) {
        throw std::invalid_argument("This function can only be used for Circular domain type.");
    }

    // 计算 y 的值
    double dx = x - circleCenterX;
    double discriminant = circleRadius * circleRadius - dx * dx;

    if (discriminant < 0) {
        throw std::invalid_argument("The given x is outside the circle.");
    }

    double sqrtDiscriminant = std::sqrt(discriminant);
    double y1 = circleCenterY + sqrtDiscriminant;
    double y2 = circleCenterY - sqrtDiscriminant;

    // 检查 y1 和 y2 是否在 [minY, maxY] 范围内
    if (y1 >= minY && y1 <= maxY) {
        return y1;
    } else if (y2 >= minY && y2 <= maxY) {
        return y2;
    } else {
        throw std::invalid_argument("No valid y found within the given range.");
    }
}

double ProblemDomain::getXOnCircle(double y, double minX, double maxX) const {
    if (domtype != DomainType::Circular) {
        throw std::invalid_argument("This function can only be used for Circular domain type.");
    }

    // 计算 x 的值
    double dy = y - circleCenterY;
    double discriminant = circleRadius * circleRadius - dy * dy;

    if (discriminant < 0) {
        throw std::invalid_argument("The given y is outside the circle.");
    }

    double sqrtDiscriminant = std::sqrt(discriminant);
    double x1 = circleCenterX + sqrtDiscriminant;
    double x2 = circleCenterX - sqrtDiscriminant;

    // 检查 x1 和 x2 是否在 [minX, maxX] 范围内
    if (x1 >= minX && x1 <= maxX) {
        return x1;
    } else if (x2 >= minX && x2 <= maxX) {
        return x2;
    } else {
        throw std::invalid_argument("No valid x found within the given range.");
    }
}

std::array<double, 2> ProblemDomain::getRayCircleIntersection(double x, double y) const {
    if (domtype != DomainType::Circular) {
        throw std::invalid_argument("This function can only be used for Circular domain type.");
    }

    // 将输入的 x, y 转换为 double 类型的坐标
    double pointX = static_cast<double>(x);
    double pointY = static_cast<double>(y);

    // 计算从圆心到点 (x, y) 的方向向量
    double dx = pointX - circleCenterX;
    double dy = pointY - circleCenterY;

    // 计算方向向量的长度
    double length = std::sqrt(dx * dx + dy * dy);

    // 如果点在圆心上，抛出异常（因为射线方向无法确定）
    if (length < 1e-10) {
        throw std::invalid_argument("The point (x, y) is at the circle center. Ray direction is undefined.");
    }

    // 归一化方向向量
    double nx = dx / length;
    double ny = dy / length;

    // 计算射线与圆的交点
    double intersectionX = circleCenterX + nx * circleRadius;
    double intersectionY = circleCenterY + ny * circleRadius;

    // 返回交点的横纵坐标
    return {intersectionX, intersectionY};
}

std::array<double, 2> ProblemDomain::getRaySquareIntersection(int i, int j) const {
    if (domtype != DomainType::Circular) {
        throw std::invalid_argument("This function can only be used for Circular domain type.");
    }

    // 检查 i 和 j 是否在有效范围内
    if (i <= 0 || i >= GridNum || j <= 0 || j >= GridNum) {
        throw std::invalid_argument("i and j must be within the range (0, GridNum).");
    }

    // 计算网格步长
    double h = 1.0 / GridNum;

    // 计算点 (i * h, j * h) 的坐标
    double pointX = i * h;
    double pointY = j * h;

    // 计算从圆心到点 (i * h, j * h) 的方向向量
    double dx = pointX - circleCenterX;
    double dy = pointY - circleCenterY;

    // 计算方向向量的长度
    double length = std::sqrt(dx * dx + dy * dy);

    // 如果点在圆心上，抛出异常（因为射线方向无法确定）
    if (length < 1e-10) {
        throw std::invalid_argument("The point (i * h, j * h) is at the circle center. Ray direction is undefined.");
    }

    // 归一化方向向量
    double nx = dx / length;
    double ny = dy / length;

    // 定义方形区域的边界
    double left = (i - 1) * h;
    double right = (i + 1) * h;
    double bottom = (j - 1) * h;
    double top = (j + 1) * h;

    // 计算射线与方形区域的交点
    // 射线方程为: (x, y) = (circleCenterX, circleCenterY) + t * (nx, ny)
    // 我们需要找到 t，使得 (x, y) 在方形区域的边界上

    // 计算与左右边界的交点
    double tLeft = (left - circleCenterX) / nx;
    double tRight = (right - circleCenterX) / nx;

    // 计算与上下边界的交点
    double tBottom = (bottom - circleCenterY) / ny; 
    double tTop = (top - circleCenterY) / ny;

    // 找到最小的正 t，即第一个与方形区域边界的交点
    double tMin = std::numeric_limits<double>::max();
    std::array<double, 2> intersection = {0.0, 0.0};

    auto checkIntersection = [&](double t) {
        if (t > 0) {
            double x = circleCenterX + t * nx;
            double y = circleCenterY + t * ny;
            if (x >= left && x <= right && y >= bottom && y <= top && !isPointInsideCircle(x, y)) {
                if (t < tMin) {
                    tMin = t;
                    intersection = {x, y};
                }
            }
        }
    };

    checkIntersection(tLeft);
    checkIntersection(tRight);
    checkIntersection(tBottom);
    checkIntersection(tTop);

    // 如果没有找到交点，抛出异常
    if (tMin == std::numeric_limits<double>::max()) {
        throw std::runtime_error("No intersection found within the square region.");
    }

    return intersection;
}

void ProblemDomain::SetValueAt00(double value) {
    valueAt00 = value;
}

double ProblemDomain::getValueAt00() const{
    return valueAt00;
}


// 导出边界点到 CSV 文件
void ProblemDomain::exportBoundaryPointsToCSV(const std::string& filename) const {
    double h = 1.0 / GridNum;
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // 遍历所有格点
    for (int j = 0; j <= GridNum; j++) {
        for (int i = 0; i <= GridNum; i++) {
            double x = i * h;
            double y = j * h;

            // 判断当前点是否在圆外
            if (!isPointInsideCircle(x, y)) {
                // 检查相邻点是否在圆内
                bool isBoundaryPoint = false;
                if (i > 0 && isPointInsideCircle((i - 1) * h, j * h)) isBoundaryPoint = true; // 左
                if (i < GridNum && isPointInsideCircle((i + 1) * h, j * h)) isBoundaryPoint = true; // 右
                if (j > 0 && isPointInsideCircle(i * h, (j - 1) * h)) isBoundaryPoint = true; // 下
                if (j < GridNum && isPointInsideCircle(i * h, (j + 1) * h)) isBoundaryPoint = true; // 上

                // 输出 1 或 0
                outFile << (isBoundaryPoint ? "1" : "0");
            } else {
                // 当前点在圆内，输出 0
                outFile << "0";
            }

            // 添加逗号（除了最后一列）
            if (i < GridNum) outFile << ",";
        }
        // 换行（除了最后一行）
        if (j < GridNum) outFile << std::endl;
    }

    outFile.close();
    std::cout << "Boundary points exported to " << filename << std::endl;
}