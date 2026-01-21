#ifndef POLYNOMIALS_AND_EXPONENTIAL_FUNCTIONS_HPP
#define POLYNOMIALS_AND_EXPONENTIAL_FUNCTIONS_HPP

#include "expression.hpp"

class VariableX : public Exp {
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class Constant : public Exp {   // Rule #1: (c*f)' = c*f'
    public:
        double value;
        bool hasFraction = false;
        long long num = 0;
        long long den = 1;
        Constant(double v);
        Constant(long long n, long long d);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};


class Power : public Exp {  // Rule #2: (x^n)' = n*x^(n-1)
    public:
        double exponent;
        bool hasFraction = false;
        long long num = 0;
        long long den = 1;
        Power(double n);
        Power(long long n, long long d);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class Exponential : public Exp { // Rule #3: (e^(a*x))' = a*e^(a*x)
    public:
        double coefficient;
        Exponential(double a);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class AddSub : public Exp {  // Rule #4: (f ± g)' = f' ± g'
    public:
        shared_ptr<Exp> left, right;
        char op;
        AddSub(shared_ptr<Exp> l, shared_ptr<Exp> r, char o);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class Multiply : public Exp { // Rule #5: (f*g)' = f'*g + f*g'
    public:
        shared_ptr<Exp> left, right;
        Multiply(shared_ptr<Exp> l, shared_ptr<Exp> r);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class Divide : public Exp { // Rule #6: (f/g)' = (f'*g - f*g')/g^2
    public:
        shared_ptr<Exp> left, right;
        Divide(shared_ptr<Exp> l, shared_ptr<Exp> r);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

#endif
