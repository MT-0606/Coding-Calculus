#ifndef CHAIN_RULE_CPP
#define CHAIN_RULE_CPP

#include "chain_rule.hpp"

#include "expression_utils.hpp"
#include "inverse_trigonometric_functions.hpp"
#include "polynomials_and_exponential_functions.hpp"
#include "trigonometric_functions.hpp"

#include <cmath>

using namespace std;

ChainRule::ChainRule(shared_ptr<Exp> f, shared_ptr<Exp> g) : outer(f), inner(g) {}
string ChainRule::toString() const {
    return "f(" + inner->toString() + ")";
}
dExp ChainRule::derivative() const {
    auto outer_deriv = outer->derivative();
    auto outer_deriv_at_g = outer_deriv->substitute(inner);
    auto inner_deriv = inner->derivative();

    return make_unique<Multiply>(
        shared_ptr<Exp>(move(outer_deriv_at_g)),
        shared_ptr<Exp>(move(inner_deriv))
    )->simplify();
}
dExp ChainRule::simplify() const {
    auto o = outer->simplify();
    auto i = inner->simplify();
    return make_unique<ChainRule>(shared_ptr<Exp>(move(o)), shared_ptr<Exp>(move(i)));
}
double ChainRule::evaluate(double x) const {
    double inner_val = inner->evaluate(x);
    return outer->evaluate(inner_val);
}
dExp ChainRule::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<ChainRule>(outer, inner->substitute(replacement))->simplify();
}

SineComposed::SineComposed(shared_ptr<Exp> a) : arg(a) {}
string SineComposed::toString() const {
    return "sin(" + arg->toString() + ")";
}
dExp SineComposed::derivative() const {
    return make_unique<Multiply>(
        make_shared<CosineComposed>(arg),
        arg->derivative()
    )->simplify();
}
dExp SineComposed::simplify() const {
    auto a = arg->simplify();
    if (dynamic_cast<VariableX*>(a.get())) {
        return make_unique<Sine>();
    }
    return make_unique<SineComposed>(shared_ptr<Exp>(move(a)));
}
double SineComposed::evaluate(double x) const {
    return sin(arg->evaluate(x));
}
dExp SineComposed::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<SineComposed>(arg->substitute(replacement))->simplify();
}

CosineComposed::CosineComposed(shared_ptr<Exp> a) : arg(a) {}
string CosineComposed::toString() const {
    return "cos(" + arg->toString() + ")";
}
dExp CosineComposed::derivative() const {
    return make_unique<Multiply>(
        make_unique<Multiply>(
            make_unique<Constant>(-1),
            make_shared<SineComposed>(arg)
        ),
        arg->derivative()
    )->simplify();
}
dExp CosineComposed::simplify() const {
    auto a = arg->simplify();
    if (dynamic_cast<VariableX*>(a.get())) {
        return make_unique<Cosine>();
    }
    return make_unique<CosineComposed>(shared_ptr<Exp>(move(a)));
}
double CosineComposed::evaluate(double x) const {
    return cos(arg->evaluate(x));
}
dExp CosineComposed::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<CosineComposed>(arg->substitute(replacement))->simplify();
}

PowerComposed::PowerComposed(shared_ptr<Exp> a, double n) : arg(a), exponent(n) {}
PowerComposed::PowerComposed(shared_ptr<Exp> a, long long n, long long d)
    : arg(a), exponent(static_cast<double>(n) / static_cast<double>(d)) {
    hasFraction = true;
    num = n;
    den = d;
    normaliseFraction(num, den);
}
string PowerComposed::toString() const {
    if (hasFraction) {
        if (den == 1) return "(" + arg->toString() + ")^" + to_string(num);
        return "(" + arg->toString() + ")^(" + formatFraction(num, den) + ")";
    }
    return "(" + arg->toString() + ")^" + formatNumber(exponent);
}
dExp PowerComposed::derivative() const {
    if (hasFraction) {
        long long n = num;
        long long d = den;
        long long n_minus = n - d;
        return make_unique<Multiply>(
            make_unique<Multiply>(
                make_unique<Constant>(n, d),
                make_shared<PowerComposed>(arg, n_minus, d)
            ),
            arg->derivative()
        )->simplify();
    }
    return make_unique<Multiply>(
        make_unique<Multiply>(
            make_unique<Constant>(exponent),
            make_shared<PowerComposed>(arg, exponent - 1)
        ),
        arg->derivative()
    )->simplify();
}
dExp PowerComposed::simplify() const {
    auto a = arg->simplify();
    if (hasFraction) {
        if (num == 0) return make_unique<Constant>(1);
        if (den == 1 && num == 1) return a;
        if (dynamic_cast<VariableX*>(a.get())) {
            return make_unique<Power>(num, den);
        }
        if (auto c = dynamic_cast<Constant*>(a.get())) {
            return make_unique<Constant>(pow(c->value, exponent));
        }
        return make_unique<PowerComposed>(shared_ptr<Exp>(move(a)), num, den);
    }
    if (exponent == 0.0) return make_unique<Constant>(1);
    if (exponent == 1.0) return a;
    if (dynamic_cast<VariableX*>(a.get())) {
        return make_unique<Power>(exponent);
    }
    if (auto c = dynamic_cast<Constant*>(a.get())) {
        return make_unique<Constant>(pow(c->value, exponent));
    }
    return make_unique<PowerComposed>(shared_ptr<Exp>(move(a)), exponent);
}
double PowerComposed::evaluate(double x) const {
    return pow(arg->evaluate(x), exponent);
}
dExp PowerComposed::substitute(const shared_ptr<Exp>& replacement) const {
    if (hasFraction) {
        return make_unique<PowerComposed>(arg->substitute(replacement), num, den)->simplify();
    }
    return make_unique<PowerComposed>(arg->substitute(replacement), exponent)->simplify();
}

ExponentialComposed::ExponentialComposed(shared_ptr<Exp> a) : arg(a) {}
string ExponentialComposed::toString() const {
    return "e^(" + arg->toString() + ")";
}
dExp ExponentialComposed::derivative() const {
    return make_unique<Multiply>(
        make_shared<ExponentialComposed>(arg),
        arg->derivative()
    )->simplify();
}
dExp ExponentialComposed::simplify() const {
    auto a = arg->simplify();
    if (dynamic_cast<VariableX*>(a.get())) {
        return make_unique<Exponential>(1);
    }
    if (auto c = dynamic_cast<Constant*>(a.get())) {
        return make_unique<Constant>(exp(c->value));
    }
    return make_unique<ExponentialComposed>(shared_ptr<Exp>(move(a)));
}
double ExponentialComposed::evaluate(double x) const {
    return exp(arg->evaluate(x));
}
dExp ExponentialComposed::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<ExponentialComposed>(arg->substitute(replacement))->simplify();
}

#endif