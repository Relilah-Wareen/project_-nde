#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include <vector>

typedef std::vector<double> Vector;

template <int Dim>
class Operator {
public:
    virtual ~Operator() = default;
    virtual Vector operator()(const Vector& v) const = 0;
};

template <int Dim>
class Injection : public Operator<Dim> {
public:
    Vector operator()(const Vector& v) const override;
};

template <int Dim>
class LinearInterpolation : public Operator<Dim> {
public:
    Vector operator()(const Vector& v) const override;
};

template <int Dim>
class FullWeighting : public Operator<Dim> {
public:
    Vector operator()(const Vector& v) const override;
};

template <int Dim>
class QuadraticInterpolation : public Operator<Dim> {
public:
    Vector operator()(const Vector& v) const override;
};


#endif // OPERATOR_HPP