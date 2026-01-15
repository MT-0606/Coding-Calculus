#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <sstream>
#include <iomanip>
using namespace std;

static inline string formatNumber(double v) {
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return (v > 0) ? "inf" : "-inf";

    double intpart;
    if (std::modf(v, &intpart) == 0.0) {
        // print as integer if the value is integral
        return std::to_string(static_cast<long long>(intpart));
    }

    // otherwise print with a reasonable precision and strip trailing zeros
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(8) << v; // adjust precision if needed
    std::string s = ss.str();
    // strip trailing zeros and possibly the trailing decimal point
    if (s.find('.') != std::string::npos) {
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
    }
    return s;
}

class Exp {    // Parent class: the expression to differentiate
    public:
        virtual ~Exp() = default;
        virtual string toString() const = 0;    // displays expressions and their derivatives
        virtual unique_ptr<Exp> derivative() const = 0; // calculates derivative of expression
        virtual unique_ptr<Exp> simplify() const = 0; // simplifies expression
        virtual double evaluate(double x) const = 0; // evaluates expression at given x
    };

using dExp = unique_ptr<Exp>;

class Constant : public Exp { // Child class #1: constant value c
    public:
        double value;
        Constant(double v);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
};

class Variable : public Exp { // Child class #2: variable x
    public:
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
};

class Power : public Exp { // Rule #1: (x^n)' = n*x^(n-1)
    public:
        double exponent;
        Power(double n);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
};

class Exponential : public Exp {   // Rule #2: (e^(a*x))' = a*e^(a*x)
    public:
        double coefficient;
        Exponential(double a);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
};

class AddSub : public Exp {    // Rule #3: (f±g)' = f'±g'
    public:
        shared_ptr<Exp> left, right;
        char op;
        AddSub(shared_ptr<Exp> l, shared_ptr<Exp> r, char o);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
};

class Multiply : public Exp {    // Rule #4: (f*g)' = f'*g+f*g'
    public:
        shared_ptr<Exp> left, right;
        Multiply(shared_ptr<Exp> l, shared_ptr<Exp> r);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
};

class Divide : public Exp {    // Rule #5: (f/g)' = (f'*g-f*g')/(g^2)
    public:
        shared_ptr<Exp> left, right;
        Divide(shared_ptr<Exp> l, shared_ptr<Exp> r);
        string toString() const override;
        dExp derivative() const override;
        dExp simplify() const override;
        double evaluate(double x) const override;
};