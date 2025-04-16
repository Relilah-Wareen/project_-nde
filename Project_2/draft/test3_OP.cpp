#include "Operator.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

// ��ӡ����
void printVector(const Vector& v, const std::string& name) {
    std::cout << name << " (" << v.size() << "): [";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << std::setw(8) << std::setprecision(4) << v[i];
        if (i < v.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

// ����FullWeighting
void testFullWeighting() {
    std::cout << "=== ����FullWeighting���� ===\n";
    FullWeighting<1> fw;
    
    // ����5����(2^2+1)
    Vector v5 = {1.0, 2.0, 4.0, 8.0, 16.0};
    Vector result = fw(v5);
    printVector(v5, "��������");
    printVector(result, "������");
    std::cout << "���ۼ���:\n";
    std::cout << "result[0] = v[0] = 1.0\n";
    std::cout << "result[1] = 0.25*v[1] + 0.5*v[2] + 0.25*v[3] = 0.25*2 + 0.5*4 + 0.25*8 = 3.0\n";
    std::cout << "result[2] = v[4] = 16.0\n\n";
    
    // ����9����(2^3+1)
    Vector v9 = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0};
    result = fw(v9);
    printVector(v9, "��������");
    printVector(result, "������");
    std::cout << "���ۼ���:\n";
    std::cout << "result[0] = v[0] = 1.0\n";
    std::cout << "result[1] = 0.25*v[1] + 0.5*v[2] + 0.25*v[3] = 3.0\n";
    std::cout << "result[2] = 0.25*v[3] + 0.5*v[4] + 0.25*v[5] = 12.0\n";
    std::cout << "result[3] = 0.25*v[5] + 0.5*v[6] + 0.25*v[7] = 48.0\n";
    std::cout << "result[4] = v[8] = 256.0\n\n";
}

// ����QuadraticInterpolation
void testQuadraticInterpolation() {
    std::cout << "=== ����QuadraticInterpolation���� ===\n";
    QuadraticInterpolation<1> quad;
    
    // ����3����(2^1+1) - Ӧ���˵����Բ�ֵ
    Vector v3 = {1.0, 2.0, 3.0};
    Vector result = quad(v3);
    printVector(v3, "��������");
    printVector(result, "������");
    std::cout << "������2ʱ���˵����Բ�ֵ:\n";
    std::cout << "result[1] = (v[0]+v[1])/2 = 1.5\n";
    std::cout << "result[3] = (v[1]+v[2])/2 = 2.5\n\n";
    
    // ����5����(2^2+1)
    Vector v5 = {1.0, 0.0, 1.0, 0.0, 1.0};
    result = quad(v5);
    printVector(v5, "��������");
    printVector(result, "������");
    std::cout << "���ۼ���:\n";
    std::cout << "result[1] = 0.375*v[0] + 0.75*v[1] - 0.125*v[2] = 0.375*1 + 0.75*0 - 0.125*1 = 0.25\n";
    std::cout << "result[3] = 9/16*(v[1]+v[2]) - 1/16*(v[0]+v[3]) = 9/16*(0+1) - 1/16*(1+0) = 0.5\n";
    std::cout << "result[5] = 9/16*(v[2]+v[3]) - 1/16*(v[1]+v[4]) = 9/16*(1+0) - 1/16*(0+1) = 0.5\n";
    std::cout << "result[7] = 0.375*v[4] + 0.75*v[3] - 0.125*v[2] = 0.375*1 + 0.75*0 - 0.125*1 = 0.25\n\n";
    
    // ����9����(2^3+1)
    Vector v9 = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    result = quad(v9);
    printVector(v9, "��������");
    printVector(result, "������");
    std::cout << "������ķ�ֵ��Χ�Ķ��β�ֵЧ��\n\n";
}

int main() {
    testFullWeighting();
    testQuadraticInterpolation();
    return 0;
}