#ifndef CHAIN_RULE_HPP
#define CHAIN_RULE_HPP

#include "expression.hpp"

class ChainRule : public Exp {
    public:
        shared_ptr<Exp> outer;
        shared_ptr<Exp> inner;
        ChainRule(shared_ptr<Exp> f, shared_ptr<Exp> g);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class SineComposed : public Exp {
    public:
        shared_ptr<Exp> arg;
        explicit SineComposed(shared_ptr<Exp> a);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class CosineComposed : public Exp {
    public:
        shared_ptr<Exp> arg;
        explicit CosineComposed(shared_ptr<Exp> a);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class PowerComposed : public Exp {
    public:
        shared_ptr<Exp> arg;
        double exponent;
        bool hasFraction = false;
        long long num = 0;
        long long den = 1;
        PowerComposed(shared_ptr<Exp> a, double n);
        PowerComposed(shared_ptr<Exp> a, long long n, long long d);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class ExponentialComposed : public Exp {
    public:
        shared_ptr<Exp> arg;
        explicit ExponentialComposed(shared_ptr<Exp> a);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

#endif
