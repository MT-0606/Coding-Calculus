#ifndef TRIGONOMETRIC_FUNCTIONS_CPP
#define TRIGONOMETRIC_FUNCTIONS_CPP

#include "trigonometric_functions.hpp"
#include "chain_rule.hpp"
#include "polynomials_and_exponential_functions.hpp"

#include <cmath>

using namespace std;

string Sine::toString() const {
    return "sin(x)";
}
dExp Sine::derivative() const {
    return make_unique<Cosine>();
}
dExp Sine::simplify() const {
    return make_unique<Sine>();
}
double Sine::evaluate(double x) const {
    return sin(x);
}

string Cosine::toString() const {
    return "cos(x)";
}
dExp Cosine::derivative() const {
    return make_unique<Multiply>(
        make_shared<Constant>(-1),
        make_shared<Sine>()
    )->simplify();
}
dExp Cosine::simplify() const {
    return make_unique<Cosine>();
}
double Cosine::evaluate(double x) const {
    return cos(x);
}

string Tangent::toString() const {
    return "tan(x)";
}
dExp Tangent::derivative() const {
    return make_unique<Multiply>(
        make_shared<Secant>(),
        make_shared<Secant>()
    )->simplify();
}
dExp Tangent::simplify() const {
    return make_unique<Tangent>();
}
double Tangent::evaluate(double x) const {
    return tan(x);
}

string Cosecant::toString() const {
    return "csc(x)";
}
dExp Cosecant::derivative() const {
    return make_unique<Multiply>(
        make_shared<Constant>(-1),
        make_shared<Multiply>(
            make_shared<Cosecant>(),
            make_shared<Cotangent>()
        )
    )->simplify();
}
dExp Cosecant::simplify() const {
    return make_unique<Cosecant>();
}
double Cosecant::evaluate(double x) const {
    return 1.0 / sin(x);
}

string Secant::toString() const {
    return "sec(x)";
}
dExp Secant::derivative() const {
    return make_unique<Multiply>(
        make_shared<Secant>(),
        make_shared<Tangent>()
    )->simplify();
}
dExp Secant::simplify() const {
    return make_unique<Secant>();
}
double Secant::evaluate(double x) const {
    return 1.0 / cos(x);
}

string Cotangent::toString() const {
    return "cot(x)";
}
dExp Cotangent::derivative() const {
    return make_unique<Multiply>(
        make_shared<Constant>(-1),
        make_unique<Divide>(
            make_shared<Constant>(1),
            make_unique<Multiply>(
                make_shared<Sine>(),
                make_shared<Sine>()
            )
        )
    )->simplify();
}
dExp Cotangent::simplify() const {
    return make_unique<Cotangent>();
}
double Cotangent::evaluate(double x) const {
    return 1.0 / tan(x);
}

dExp Sine::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<SineComposed>(replacement)->simplify();
}
dExp Cosine::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<CosineComposed>(replacement)->simplify();
}
dExp Tangent::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<Tangent>();
}
dExp Cosecant::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<Divide>(
        make_shared<Constant>(1),
        make_shared<SineComposed>(replacement)
    )->simplify();
}
dExp Secant::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<Divide>(
        make_shared<Constant>(1),
        make_shared<CosineComposed>(replacement)
    )->simplify();
}
dExp Cotangent::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<Divide>(
        make_shared<CosineComposed>(replacement),
        make_shared<SineComposed>(replacement)
    )->simplify();
}

#endif