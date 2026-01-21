#ifndef INVERSE_TRIGONOMETRIC_FUNCTIONS_CPP
#define INVERSE_TRIGONOMETRIC_FUNCTIONS_CPP

#include "inverse_trigonometric_functions.hpp"
#include "polynomials_and_exponential_functions.hpp"

#include <cmath>

using namespace std;

Sqrt::Sqrt(shared_ptr<Exp> a) : arg(a) {}
string Sqrt::toString() const {
    return "sqrt(" + arg->toString() + ")";
}
dExp Sqrt::derivative() const {
    return make_unique<Divide>(
        arg->derivative(),
        make_unique<Multiply>(
            make_shared<Constant>(2),
            make_shared<Sqrt>(arg)
        )
    )->simplify();
}
dExp Sqrt::simplify() const {
    auto a = arg->simplify();
    if (auto c = dynamic_cast<Constant*>(a.get())) {
        if (c->value >= 0) return make_unique<Constant>(sqrt(c->value));
    }
    return make_unique<Sqrt>(shared_ptr<Exp>(move(a)));
}
double Sqrt::evaluate(double x) const {
    return sqrt(arg->evaluate(x));
}

string ArcSine::toString() const {
    return "arcsin(x)";
}
dExp ArcSine::derivative() const {
    return make_unique<Divide>(
        make_unique<Constant>(1),
        make_unique<Sqrt>(
            make_shared<AddSub>(
                make_shared<Constant>(1),
                make_shared<Power>(2),
                '-'
            )
        )
    )->simplify();
}
dExp ArcSine::simplify() const {
    return make_unique<ArcSine>();
}
double ArcSine::evaluate(double x) const {
    return asin(x);
}

string ArcCosine::toString() const {
    return "arccos(x)";
}
dExp ArcCosine::derivative() const {
    return make_unique<Divide>(
        make_unique<Constant>(-1),
        make_unique<Sqrt>(
            make_shared<AddSub>(
                make_shared<Constant>(1),
                make_shared<Power>(2),
                '-'
            )
        )
    )->simplify();
}
dExp ArcCosine::simplify() const {
    return make_unique<ArcCosine>();
}
double ArcCosine::evaluate(double x) const {
    return acos(x);
}

string ArcTangent::toString() const {
    return "arctan(x)";
}
dExp ArcTangent::derivative() const {
    return make_unique<Divide>(
        make_unique<Constant>(1),
        make_unique<AddSub>(
            make_shared<Constant>(1),
            make_shared<Power>(2),
            '+'
        )
    )->simplify();
}
dExp ArcTangent::simplify() const {
    return make_unique<ArcTangent>();
}
double ArcTangent::evaluate(double x) const {
    return atan(x);
}

string ArcCosecant::toString() const {
    return "arccsc(x)";
}
dExp ArcCosecant::derivative() const {
    auto absx = make_shared<Sqrt>(make_shared<Power>(2));
    auto root = make_shared<Sqrt>(
        make_shared<AddSub>(
            make_shared<Power>(2),
            make_shared<Constant>(1),
            '-'
        )
    );
    return make_unique<Divide>(
        make_unique<Constant>(-1),
        make_unique<Multiply>(absx, root)
    )->simplify();
}
dExp ArcCosecant::simplify() const {
    return make_unique<ArcCosecant>();
}
double ArcCosecant::evaluate(double x) const {
    return asin(1.0 / x);
}

string ArcSecant::toString() const {
    return "arcsec(x)";
}
dExp ArcSecant::derivative() const {
    auto absx = make_shared<Sqrt>(make_shared<Power>(2));
    auto root = make_shared<Sqrt>(
        make_shared<AddSub>(
            make_shared<Power>(2),
            make_shared<Constant>(1),
            '-'
        )
    );
    return make_unique<Divide>(
        make_unique<Constant>(1),
        make_unique<Multiply>(absx, root)
    )->simplify();
}
dExp ArcSecant::simplify() const {
    return make_unique<ArcSecant>();
}
double ArcSecant::evaluate(double x) const {
    return acos(1.0 / x);
}

string ArcCotangent::toString() const {
    return "arccot(x)";
}
dExp ArcCotangent::derivative() const {
    return make_unique<Divide>(
        make_unique<Constant>(-1),
        make_unique<AddSub>(
            make_shared<Constant>(1),
            make_shared<Power>(2),
            '+'
        )
    )->simplify();
}
dExp ArcCotangent::simplify() const {
    return make_unique<ArcCotangent>();
}
double ArcCotangent::evaluate(double x) const {
    return atan(1.0 / x);
}

dExp Sqrt::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<Sqrt>(arg->substitute(replacement))->simplify();
}
dExp ArcSine::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<ArcSine>();
}
dExp ArcCosine::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<ArcCosine>();
}
dExp ArcTangent::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<ArcTangent>();
}
dExp ArcCosecant::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<ArcCosecant>();
}
dExp ArcSecant::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<ArcSecant>();
}
dExp ArcCotangent::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<ArcCotangent>();
}

#endif