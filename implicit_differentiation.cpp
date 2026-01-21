#ifndef IMPLICIT_DIFFERENTIATION_CPP
#define IMPLICIT_DIFFERENTIATION_CPP

#include "implicit_differentiation.hpp"

#include "chain_rule.hpp"
#include "polynomials_and_exponential_functions.hpp"
#include "trigonometric_functions.hpp"
#include "inverse_trigonometric_functions.hpp"

#include <cmath>

using namespace std;

string VariableY::toString() const {
    return "y"; 
}
dExp VariableY::derivative() const {
    return make_unique<DerivativeY>();
}
dExp VariableY::simplify() const {
    return make_unique<VariableY>();
}
double VariableY::evaluate(double x) const {
    return NAN;
}
dExp VariableY::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<VariableY>();
}

string DerivativeY::toString() const {
    return "y'";
}
dExp DerivativeY::derivative() const {
    return make_unique<DerivativeY>();
}
dExp DerivativeY::simplify() const {
    return make_unique<DerivativeY>();
}
double DerivativeY::evaluate(double x) const { 
    return NAN;
}
dExp DerivativeY::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<DerivativeY>();
}

ImplicitEquation::ImplicitEquation(shared_ptr<Exp> l, shared_ptr<Exp> r) : left(l), right(r) {}
string ImplicitEquation::toString() const {
    return left->toString() + " = " + right->toString();
}

static bool containsYPrime(const Exp* expr) {
    if (dynamic_cast<const DerivativeY*>(expr)) return true;
    if (auto add = dynamic_cast<const AddSub*>(expr)) {
        return containsYPrime(add->left.get()) || containsYPrime(add->right.get());
    }
    if (auto mul = dynamic_cast<const Multiply*>(expr)) {
        return containsYPrime(mul->left.get()) || containsYPrime(mul->right.get());
    }
    if (auto div = dynamic_cast<const Divide*>(expr)) {
        return containsYPrime(div->left.get()) || containsYPrime(div->right.get());
    }
    return false;
}

static shared_ptr<Exp> asShared(dExp expr) {
    return shared_ptr<Exp>(move(expr));
}

static dExp makeZero() { return make_unique<Constant>(0); }
static dExp makeOne() { return make_unique<Constant>(1); }

static dExp addExpr(dExp a, dExp b, char op) {
    return make_unique<AddSub>(asShared(move(a)), asShared(move(b)), op)->simplify();
}
static dExp mulExpr(dExp a, dExp b) {
    return make_unique<Multiply>(asShared(move(a)), asShared(move(b)))->simplify();
}
static dExp divExpr(dExp a, dExp b) {
    return make_unique<Divide>(asShared(move(a)), asShared(move(b)))->simplify();
}

static bool isZeroConst(const dExp& expr) {
    if (auto c = dynamic_cast<Constant*>(expr.get())) {
        return c->value == 0.0;
    }
    return false;
}

// Split expr into a*Y' + b where a,b are expressions without Y'
static bool splitLinearYPrime(const shared_ptr<Exp>& expr, dExp& coeff, dExp& rest) {
    if (dynamic_cast<DerivativeY*>(expr.get())) {
        coeff = makeOne();
        rest = makeZero();
        return true;
    }
    if (!containsYPrime(expr.get())) {
        coeff = makeZero();
        rest = expr->simplify();
        return true;
    }
    if (auto add = dynamic_cast<AddSub*>(expr.get())) {
        dExp lc, lr, rc, rr;
        if (!splitLinearYPrime(add->left, lc, lr)) return false;
        if (!splitLinearYPrime(add->right, rc, rr)) return false;
        coeff = addExpr(move(lc), move(rc), add->op);
        rest = addExpr(move(lr), move(rr), add->op);
        return true;
    }
    if (auto mul = dynamic_cast<Multiply*>(expr.get())) {
        bool lHas = containsYPrime(mul->left.get());
        bool rHas = containsYPrime(mul->right.get());
        if (lHas && rHas) return false;

        if (lHas) {
            dExp lc, lr;
            if (!splitLinearYPrime(mul->left, lc, lr)) return false;
            if (!isZeroConst(lr)) return false;
            coeff = mulExpr(move(lc), mul->right->simplify());
            rest = makeZero();
            return true;
        }
        dExp rc, rr;
        if (!splitLinearYPrime(mul->right, rc, rr)) return false;
        if (!isZeroConst(rr)) return false;
        coeff = mulExpr(move(rc), mul->left->simplify());
        rest = makeZero();
        return true;
    }
    if (auto div = dynamic_cast<Divide*>(expr.get())) {
        if (containsYPrime(div->right.get())) return false;
        dExp nc, nr;
        if (!splitLinearYPrime(div->left, nc, nr)) return false;
        coeff = divExpr(move(nc), div->right->simplify());
        rest = divExpr(move(nr), div->right->simplify());
        return true;
    }

    return false;
}

dExp ImplicitEquation::derivative() const {
    auto dl = left->derivative();
    auto dr = right->derivative();
    auto diff = make_unique<AddSub>(asShared(move(dl)), asShared(move(dr)), '-')->simplify();

    dExp coeff;
    dExp rest;
    if (!splitLinearYPrime(asShared(move(diff)), coeff, rest)) {
        return make_unique<Constant>(NAN);
    }

    // a*Y' + b = 0 -> Y' = -b/a
    dExp negRest = mulExpr(make_unique<Constant>(-1), move(rest));
    return divExpr(move(negRest), move(coeff));
}

#endif