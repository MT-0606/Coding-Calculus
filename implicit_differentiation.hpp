#ifndef IMPLICIT_DIFFERENTIATION_HPP
#define IMPLICIT_DIFFERENTIATION_HPP

#include "expression.hpp"

class VariableY : public Exp {
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class DerivativeY : public Exp {
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class ImplicitEquation {
    public:
        shared_ptr<Exp> left;
        shared_ptr<Exp> right;
        ImplicitEquation(shared_ptr<Exp> l, shared_ptr<Exp> r);
        string toString() const;
        dExp derivative() const;
};

#endif
