#ifndef INVERSE_TRIGONOMETRIC_FUNCTIONS_HPP
#define INVERSE_TRIGONOMETRIC_FUNCTIONS_HPP

#include "expression.hpp"

class Sqrt : public Exp {
    public:
        shared_ptr<Exp> arg;
        explicit Sqrt(shared_ptr<Exp> a);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class ArcSine : public Exp {
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class ArcCosine : public Exp {
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class ArcTangent : public Exp {
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class ArcCosecant : public Exp {
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class ArcSecant : public Exp {
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

class ArcCotangent : public Exp {
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
        dExp substitute(const shared_ptr<Exp>& replacement) const override;
};

#endif
