#include "BVPSolver.hpp"

BVPSolver::BVPSolver(ProblemDomain pd, double (*function)(double, double))
    : pd(pd), GridNum(pd.getGridNum()), function(function) {
    int n = pd.getGridNum();
    F.resize((n + 1) * (n + 1), 0.0);
    A.resize((n + 1) * (n + 1), std::vector<double>((n + 1) * (n + 1), 0.0));
    h = 1.0 / GridNum;
}

void BVPSolver::ShowMatrix() {
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            std::cout << A[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

void BVPSolver::ShowVectorF() {
    for (int i = 0; i < F.size(); i++) {
        std::cout << F[i] << "\t";
    }
    std::cout << std::endl;
}

void BVPSolver::ShowMatrix(const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            outFile << A[i][j];
            if (j < A[0].size() - 1) {
                outFile << ",";
            }
        }
        outFile << std::endl;
    }

    outFile.close();
    std::cout << "Matrix written to " << filename << " in CSV format." << std::endl;
}

void BVPSolver::ShowVectorF(const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    for (int i = 0; i < F.size(); i++) {
        outFile << F[i];
        if (i < F.size() - 1) {
            outFile << "," << std::endl;
        }
    }

    outFile.close();
    std::cout << "Vector written to " << filename << " in CSV format." << std::endl;
}

void BVPSolver::BuildLineraSystem() {
    // 检查是否需要设置额外的约束条件
    if ((pd.getDomainType() == DomainType::Default && pd.getsquareBoundaryType() == BCType::Neumann) || 
    (pd.getDomainType() == DomainType::Circular && pd.getsquareBoundaryType() != BCType::Dirichlet && pd.getcircleBoundaryType() == BCType::Neumann)) {
        // 检查是否设置了 (0,0) 处的值
        if (pd.getValueAt00() == std::numeric_limits<double>::max()) {
            throw std::runtime_error("Please set the value at (0,0) using SetValueAt00().");
        }
    }
    
    int idx = 0;
    for (int j = 0; j <= GridNum; j++) {
        for (int i = 0; i <= GridNum; i++) {
            idx = i + j * (GridNum + 1);
            if (i == 0 || j == 0 || i == GridNum || j == GridNum) {
                if (pd.getsquareBoundaryType() == BCType::Dirichlet) {
                    A[idx][i + j * (GridNum + 1)] = 1;
                    if(i==0) F[idx] = pd.getSquareBoundaryValue(0, i * h, j * h);
                    else if(i==GridNum) F[idx] = pd.getSquareBoundaryValue(1, i * h, j * h);
                    else if(j==0) F[idx] = pd.getSquareBoundaryValue(2, i * h, j * h);
                    else F[idx] = pd.getSquareBoundaryValue(3, i * h, j * h);
                } 
                else if (pd.getsquareBoundaryType() == BCType::Neumann) {
                    A[idx][i + j * (GridNum + 1)] = 4;
                    if (i == 0 && j == 0) {
                        A[idx][i + (j + 1) * (GridNum + 1)] = -2;
                        A[idx][i + 1 + j * (GridNum + 1)] = -2;
                        F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(0, i * h, j * h) + 2 * h * pd.getSquareBoundaryValue(2, i * h, j * h);
                    } else if (i == 0 && j == GridNum) {
                        A[idx][i + (j - 1) * (GridNum + 1)] = -2;
                        A[idx][i + 1 + j * (GridNum + 1)] = -2;
                        F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(0, i * h, j * h) + 2 * h * pd.getSquareBoundaryValue(3, i * h, j * h);
                    } else if (i == GridNum && j == 0) {
                        A[idx][i - 1 + j * (GridNum + 1)] = -2;
                        A[idx][i + (j + 1) * (GridNum + 1)] = -2;
                        F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(1, i * h, j * h) + 2 * h * pd.getSquareBoundaryValue(2, i * h, j * h);
                    } else if (i == GridNum && j == GridNum) {
                        A[idx][i - 1 + j * (GridNum + 1)] = -2;
                        A[idx][i + (j - 1) * (GridNum + 1)] = -2;
                        F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(1, i * h, j * h) + 2 * h * pd.getSquareBoundaryValue(3, i * h, j * h);
                    } else if (j == 0) {
                        A[idx][i + (j + 1) * (GridNum + 1)] = -2;
                        A[idx][i - 1 + j * (GridNum + 1)] = -1;
                        A[idx][i + 1 + j * (GridNum + 1)] = -1;
                        F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(2, i * h, j * h);
                    } else if (i == 0) {
                        A[idx][i + 1 + j * (GridNum + 1)] = -2;
                        A[idx][i + (j + 1) * (GridNum + 1)] = -1;
                        A[idx][i + (j - 1) * (GridNum + 1)] = -1;
                        F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(0, i * h, j * h);
                    } else if (i == GridNum) {
                        A[idx][i - 1 + j * (GridNum + 1)] = -2;
                        A[idx][i + (j - 1) * (GridNum + 1)] = -1;
                        A[idx][i + (j + 1) * (GridNum + 1)] = -1;
                        F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(1, i * h, j * h);
                    } else if (j == GridNum) {
                        A[idx][i + (j - 1) * (GridNum + 1)] = -2;
                        A[idx][i - 1 + j * (GridNum + 1)] = -1;
                        A[idx][i + 1 + j * (GridNum + 1)] = -1;
                        F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(3, i * h, j * h);
                    } 
                }
                else if (pd.getsquareBoundaryType() == BCType::Mixed) {
                    std::array<BCType, 4> boundtp=pd.getEdgeTypes(); // 左右下上
                    // A[idx][i + j * (GridNum + 1)] = 4;            // 0 1 2 3
                    if (i == 0 && j == 0) {
                        if(boundtp[0]==BCType::Dirichlet || boundtp[2]==BCType::Dirichlet){
                            int index=0;
                            if(boundtp[0]==BCType::Dirichlet)index=0;
                            else index=2;
                            A[idx][i + j * (GridNum + 1)] = 1;
                            F[idx] = pd.getSquareBoundaryValue(index, i * h, j * h);
                        }
                        else {
                            A[idx][i + j * (GridNum + 1)] = 4;
                            A[idx][i + (j + 1) * (GridNum + 1)] = -2;
                            A[idx][i + 1 + j * (GridNum + 1)] = -2;
                            F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(0, i * h, j * h) + 2 * h * pd.getSquareBoundaryValue(2, i * h, j * h);    
                        }
                    } 
                    else if (i == 0 && j == GridNum) {
                        if(boundtp[0]==BCType::Dirichlet || boundtp[3]==BCType::Dirichlet){
                            int index=0;
                            if(boundtp[0]==BCType::Dirichlet)index=0;
                            else index=3;
                            A[idx][i + j * (GridNum + 1)] = 1;
                            F[idx] = pd.getSquareBoundaryValue(index, i * h, j * h);
                        }
                        else {
                            A[idx][i + j * (GridNum + 1)] = 4;
                            A[idx][i + (j - 1) * (GridNum + 1)] = -2;
                            A[idx][i + 1 + j * (GridNum + 1)] = -2;
                            F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(0, i * h, j * h) + 2 * h * pd.getSquareBoundaryValue(3, i * h, j * h);    
                        }
                    } 
                    else if (i == GridNum && j == 0) {
                        if(boundtp[2]==BCType::Dirichlet || boundtp[1]==BCType::Dirichlet){
                            int index=0;
                            if(boundtp[2]==BCType::Dirichlet)index=2;
                            else index=1;
                            A[idx][i + j * (GridNum + 1)] = 1;
                            F[idx] = pd.getSquareBoundaryValue(index, i * h, j * h);    
                        }
                        else {
                            A[idx][i + j * (GridNum + 1)] = 4;
                            A[idx][i - 1 + j * (GridNum + 1)] = -2;
                            A[idx][i + (j + 1) * (GridNum + 1)] = -2;
                            F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(1, i * h, j * h) + 2 * h * pd.getSquareBoundaryValue(2, i * h, j * h);
                        }
                    } 
                    else if (i == GridNum && j == GridNum) {
                        if(boundtp[1]==BCType::Dirichlet || boundtp[3]==BCType::Dirichlet){
                            int index=0;
                            if(boundtp[1]==BCType::Dirichlet)index=1;
                            else index=3;
                            A[idx][i + j * (GridNum + 1)] = 1;
                            F[idx] = pd.getSquareBoundaryValue(index, i * h, j * h);        
                        }
                        else {
                            A[idx][i + j * (GridNum + 1)] = 4;
                            A[idx][i - 1 + j * (GridNum + 1)] = -2;
                            A[idx][i + (j - 1) * (GridNum + 1)] = -2;
                            F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(1, i * h, j * h) + 2 * h * pd.getSquareBoundaryValue(3, i * h, j * h);    
                        }
                    } 
                    else if (j == 0) {
                        if(boundtp[2]==BCType::Dirichlet){
                            A[idx][i + j * (GridNum + 1)] = 1;
                            F[idx] = pd.getSquareBoundaryValue(2, i * h, j * h);        
                        }
                        else {
                            A[idx][i + j * (GridNum + 1)] = 4;
                            A[idx][i + (j + 1) * (GridNum + 1)] = -2;
                            A[idx][i - 1 + j * (GridNum + 1)] = -1;
                            A[idx][i + 1 + j * (GridNum + 1)] = -1;
                            F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(2, i * h, j * h);    
                        }
                    } 
                    else if (i == 0) {
                        if(boundtp[0]==BCType::Dirichlet){
                            A[idx][i + j * (GridNum + 1)] = 1;
                            F[idx] = pd.getSquareBoundaryValue(0, i * h, j * h);            
                        }
                        else {
                            A[idx][i + j * (GridNum + 1)] = 4;
                            A[idx][i + 1 + j * (GridNum + 1)] = -2;
                            A[idx][i + (j + 1) * (GridNum + 1)] = -1;
                            A[idx][i + (j - 1) * (GridNum + 1)] = -1;
                            F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(0, i * h, j * h);    
                        }
                    } 
                    else if (i == GridNum) {
                        if(boundtp[1]==BCType::Dirichlet){
                            A[idx][i + j * (GridNum + 1)] = 1;
                            F[idx] = pd.getSquareBoundaryValue(1, i * h, j * h);                
                        }
                        else {
                            A[idx][i + j * (GridNum + 1)] = 4;
                            A[idx][i - 1 + j * (GridNum + 1)] = -2;
                            A[idx][i + (j - 1) * (GridNum + 1)] = -1;
                            A[idx][i + (j + 1) * (GridNum + 1)] = -1;
                            F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(1, i * h, j * h);    
                        }
                    } 
                    else if (j == GridNum) {
                        if(boundtp[2]==BCType::Dirichlet){
                            A[idx][i + j * (GridNum + 1)] = 1;
                            F[idx] = pd.getSquareBoundaryValue(3, i * h, j * h);                    
                        }
                        else {
                            A[idx][i + j * (GridNum + 1)] = 4;
                            A[idx][i + (j - 1) * (GridNum + 1)] = -2;
                            A[idx][i - 1 + j * (GridNum + 1)] = -1;
                            A[idx][i + 1 + j * (GridNum + 1)] = -1;
                            F[idx] = function(i * h, j * h) * h * h + 2 * h * pd.getSquareBoundaryValue(3, i * h, j * h);
                        }
                    } 
                } 
                else {
                    throw std::runtime_error("Unexpected boundary type.");
                }
            } 
            else {
                A[idx][i + j * (GridNum + 1)] = 4;
                A[idx][i - 1 + j * (GridNum + 1)] = -1;
                A[idx][i + 1 + j * (GridNum + 1)] = -1;
                A[idx][i + (j - 1) * (GridNum + 1)] = -1;
                A[idx][i + (j + 1) * (GridNum + 1)] = -1;
                F[idx] = function(i * h, j * h) * h * h;
            }
            // idx++;
        }
    }

    if(pd.getDomainType()==DomainType::Circular){
        idx = 0;
        // 这里我们要求圆在内部，不判断边界.
        for(int j = 1; j < GridNum; j++){
            for(int i = 1; i < GridNum; i++){
                idx = i + j * (GridNum + 1);
                if(pd.getcircleBoundaryType()==BCType::Dirichlet){
                    if(pd.isPointInsideCircle(i*h,j*h)){
                        for(int k = 0; k <= GridNum; k++){
                            A[idx][k]=0;
                        } 
                        A[idx][idx]=100;  // 刻意设置的大一点，本质是给圆内的点置零，对其他元素无影响
                        F[idx]=0;
                    }
                    else {
                        if(pd.isPointInsideCircle((i-1)*h, j*h)){
                            double x = pd.getXOnCircle(j*h, (i-1)*h, i*h);
                            double dist = i*h-x;
                            double theta = dist / h;
                            A[idx][i - 1 + j * (GridNum + 1)] = 0;
                            A[idx][idx] += (2.0/theta)-2;
                            A[idx][i + 1 + j * (GridNum + 1)] -= (1.0-theta)/(1.0+theta);
                            F[idx] += 2.0/theta/(1+theta)*pd.getCircleBoundaryValue(x,j*h);
                        }
                        if(pd.isPointInsideCircle((i+1)*h, j*h)){
                            double x = pd.getXOnCircle(j*h, i*h, (i+1)*h);
                            double dist = (i+1)*h-x;
                            double theta = dist / h;
                            A[idx][i + 1 + j * (GridNum + 1)] = 0;
                            A[idx][idx] += (2.0/theta)-2;
                            A[idx][i - 1 + j * (GridNum + 1)] -= (1.0-theta)/(1.0+theta);
                            F[idx] += 2.0/theta/(1+theta)*pd.getCircleBoundaryValue(x,j*h);
                        }
                        if(pd.isPointInsideCircle(i*h, (j-1)*h)){
                            double y = pd.getYOnCircle(i*h, (j-1)*h, j*h);
                            double dist = y-j*h;
                            double theta = dist / h;
                            A[idx][i + (j - 1) * (GridNum + 1)] = 0;
                            A[idx][idx] += (2.0/theta)-2;
                            A[idx][i + (j + 1) * (GridNum + 1)] -= (1.0-theta)/(1.0+theta);
                            F[idx] += 2.0/theta/(1+theta)*pd.getCircleBoundaryValue(i*h,y);
                        }
                        if(pd.isPointInsideCircle(i*h, (j+1)*h)){
                            double y = pd.getYOnCircle(i*h, j*h, (j+1)*h);
                            double dist = (j+1)*h-y;
                            double theta = dist / h;
                            A[idx][i + (j + 1) * (GridNum + 1)] = 0;
                            A[idx][idx] += (2.0/theta)-2;
                            A[idx][i + (j - 1) * (GridNum + 1)] -= (1.0-theta)/(1.0+theta);
                            F[idx] += 2.0/theta/(1+theta)*pd.getCircleBoundaryValue(i*h,y);
                        }
                    }
                }
                else if(pd.getcircleBoundaryType()==BCType::Neumann){
                    // 这里我们要求圆在内部，不判断边界.
                    if(pd.isPointInsideCircle(i*h,j*h)){
                        for(int k = 0; k <= GridNum; k++){
                            A[idx][k]=0;
                        } 
                        A[idx][idx]=1919810;  // 刻意设置的大一点，本质是给圆内的点置零，对其他元素无影响
                        F[idx]=0;
                    }
                    else {
                        if (
                            pd.isPointInsideCircle((i - 1) * h, j * h) || // 左
                            pd.isPointInsideCircle((i + 1) * h, j * h) || // 右
                            pd.isPointInsideCircle(i * h, (j - 1) * h) || // 下
                            pd.isPointInsideCircle(i * h, (j + 1) * h)    // 上
                        ) 
                        {
                            // 上下左右至少有一个点在圆内部,此处采用一阶精度的插值,作圆心到该点延长线,对该点,圆上交点和"田"内交点插值.
                            std::array<double, 2> Q=pd.getRayCircleIntersection(i*h,j*h); // 与圆的交点.
                            std::array<double, 2> T=pd.getRaySquareIntersection(i, j); // 与网格交点.
                            double QP=std::sqrt((Q[0]-i*h)*(Q[0]-i*h)+(Q[1]-j*h)*(Q[1]-j*h));
                            double PT=std::sqrt((T[0]-i*h)*(T[0]-i*h)+(T[1]-j*h)*(T[1]-j*h));
                            double RT=0,TS=0;
                            int ri=0,rj=0,si=0,sj=0;
                            if(T[1]==(j-1)*h || T[1]==(j+1)*h){
                                rj=T[1];
                                sj=T[1];
                                if(T[0]>=(i-1)*h && T[0]<=i*h) {
                                    ri=i-1;
                                    si=i;
                                }
                                else {
                                    ri=i;
                                    si=i+1;
                                }
                                RT=T[0]-ri*h;
                                TS=si*h-T[0];
                            }
                            else if(T[0]==(i-1)*h || T[0]==(i+1)*h){
                                si=T[0];
                                ri=T[0];
                                if(T[1]>=(j-1)*h && T[1]<=j*h){
                                    rj=j-1;
                                    sj=j;
                                }
                                else {
                                    rj=j;
                                    sj=j+1;
                                }
                                RT=T[1]-rj*h;
                                TS=sj*h-T[1];
                            }
                            else {
                                throw std::runtime_error("Unexpected error when build circular Neumann bound.");
                            }
                            
                            for(int k=0;k<=GridNum;k++){
                                A[idx][k]=0;
                            }
                            //idx = i + j * (GridNum + 1);
                            A[idx][idx]=(-1.0)/PT;
                            A[idx][ri+rj*(GridNum+1)] = TS/PT/h;
                            A[idx][si+sj*(GridNum+1)] = RT/PT/h;
                            F[idx] = pd.getCircleBoundaryValue(Q[0],Q[1]);
                        }
                        else {} // do nothing.
                    }
                }
                else {} // do nothing.
            }
        }
    }
    else {}  // do nothing.

    // 乘一个大的数用于改善条件数
    if ((pd.getDomainType() == DomainType::Default && pd.getsquareBoundaryType() == BCType::Neumann) || 
(pd.getDomainType() == DomainType::Circular && pd.getsquareBoundaryType() != BCType::Dirichlet && pd.getcircleBoundaryType() == BCType::Neumann)) {
        A[0][0] = A[0][0] + 100;
        F[0] = F[0] + 100*pd.getValueAt00();
    }

}

std::vector<double> BVPSolver::solve() {
    int n = A.size();
    int nrhs = 1;
    std::vector<int> ipiv(n);

    std::vector<double> A_flat;
    for (const auto& row : A) {
        A_flat.insert(A_flat.end(), row.begin(), row.end());
    }

    int info = LAPACKE_dgesv(
        LAPACK_ROW_MAJOR,
        n,
        nrhs,
        A_flat.data(),
        n,
        ipiv.data(),
        F.data(),
        nrhs
    );

    if (info == 0) {
    } else {
        std::cerr << "Error: Failed to solve the linear system. Info = " << info << std::endl;
        return F;
    }
    std::cout << "FINISH SOLVE." << std::endl;
    return F;
}


int BVPSolver::getGridNum() const{
    return GridNum;
}