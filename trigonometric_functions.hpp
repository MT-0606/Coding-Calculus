#ifndef TRIGONOMETRIC_FUNCTIONS_HPP
#define TRIGONOMETRIC_FUNCTIONS_HPP

#include "expression.hpp"

class Sine : public Exp { // (sin(x))' = cos(x)
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class Cosine : public Exp { // (cos(x))' = -sin(x)
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class Tangent : public Exp { // (tan(x))' = sec^2(x)
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class Cosecant : public Exp { // (csc(x))' = -csc(x)*cot(x)
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class Secant : public Exp { // (sec(x))' = sec(x)*tan(x)
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class Cotangent : public Exp { // (cot(x))' = -(csc(x))^2
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

#endif
