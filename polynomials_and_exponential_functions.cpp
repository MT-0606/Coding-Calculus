#ifndef POLYNOMIALS_AND_EXPONENTIAL_FUNCTIONS_CPP
#define POLYNOMIALS_AND_EXPONENTIAL_FUNCTIONS_CPP

#include "polynomials_and_exponential_functions.hpp"

#include "chain_rule.hpp"
#include "expression_utils.hpp"
#include "inverse_trigonometric_functions.hpp"
#include "trigonometric_functions.hpp"

#include <map>
#include <utility>
#include <vector>

using namespace std;

static inline shared_ptr<Exp> toShared(dExp expr) {
    return shared_ptr<Exp>(move(expr));
}
static inline Constant* asConst(const shared_ptr<Exp>& expr) {
    return dynamic_cast<Constant*>(expr.get());
}
static inline bool getRational(const Constant* c, long long& n, long long& d) {
    if (c->hasFraction) {
        n = c->num;
        d = c->den;
        normaliseFraction(n, d);
        return true;
    }
    if (isIntegerDouble(c->value)) {
        n = llround(c->value);
        d = 1;
        return true;
    }
    return false;
}
static inline dExp makeRationalConst(long long n, long long d) {
    normaliseFraction(n, d);
    return make_unique<Constant>(n, d);
}
static inline bool isConstValue(const Exp* expr, double v) {
    if (auto c = dynamic_cast<const Constant*>(expr)) {
        return fabs(c->value - v) < 1e-12;
    }
    return false;
}
static inline bool isTangentExpr(const Exp* expr) {
    return dynamic_cast<const Tangent*>(expr) != nullptr;
}
static inline bool isSecantExpr(const Exp* expr) {
    return dynamic_cast<const Secant*>(expr) != nullptr;
}
static inline bool isAddOnePlusTangent(const Exp* expr) {
    auto add = dynamic_cast<const AddSub*>(expr);
    if (!add || add->op != '+') return false;
    if (isConstValue(add->left.get(), 1.0) && isTangentExpr(add->right.get())) return true;
    if (isConstValue(add->right.get(), 1.0) && isTangentExpr(add->left.get())) return true;
    return false;
}
static inline bool isTanTimesOnePlusTan(const Exp* expr) {
    auto mul = dynamic_cast<const Multiply*>(expr);
    if (!mul) return false;
    return (isTangentExpr(mul->left.get()) && isAddOnePlusTangent(mul->right.get())) ||
           (isTangentExpr(mul->right.get()) && isAddOnePlusTangent(mul->left.get()));
}
static inline bool isSecSquaredExpr(const Exp* expr) {
    auto mul = dynamic_cast<const Multiply*>(expr);
    if (!mul) return false;
    return isSecantExpr(mul->left.get()) && isSecantExpr(mul->right.get());
}
static inline bool isAddSubOrDivideExpr(const Exp* expr) {
    return dynamic_cast<const AddSub*>(expr) || dynamic_cast<const Divide*>(expr);
}
static inline bool isPowerLikeExpr(const Exp* expr) {
    return dynamic_cast<const Power*>(expr) || dynamic_cast<const PowerComposed*>(expr);
}
static inline bool isTrigLikeExpr(const Exp* expr) {
    return dynamic_cast<const Sine*>(expr) ||
           dynamic_cast<const Cosine*>(expr) ||
           dynamic_cast<const Tangent*>(expr) ||
           dynamic_cast<const Cosecant*>(expr) ||
           dynamic_cast<const Secant*>(expr) ||
           dynamic_cast<const Cotangent*>(expr) ||
           dynamic_cast<const SineComposed*>(expr) ||
           dynamic_cast<const CosineComposed*>(expr) ||
           dynamic_cast<const ArcSine*>(expr) ||
           dynamic_cast<const ArcCosine*>(expr) ||
           dynamic_cast<const ArcTangent*>(expr) ||
           dynamic_cast<const ArcCosecant*>(expr) ||
           dynamic_cast<const ArcSecant*>(expr) ||
           dynamic_cast<const ArcCotangent*>(expr);
}
static inline bool isExponentialLikeExpr(const Exp* expr) {
    return dynamic_cast<const Exponential*>(expr) ||
           dynamic_cast<const ExponentialComposed*>(expr);
}

static void collectFactors(const shared_ptr<Exp>& expr, vector<shared_ptr<Exp>>& out) {
    if (auto mul = dynamic_cast<Multiply*>(expr.get())) {
        collectFactors(mul->left, out);
        collectFactors(mul->right, out);
        return;
    }
    out.push_back(expr);
}

static shared_ptr<Exp> buildProduct(const vector<shared_ptr<Exp>>& factors) {
    if (factors.empty()) return make_shared<Constant>(1.0);
    shared_ptr<Exp> acc = factors[0];
    for (size_t i = 1; i < factors.size(); ++i) {
        acc = make_shared<Multiply>(acc, factors[i]);
    }
    return acc;
}
static dExp buildProductUnique(const vector<shared_ptr<Exp>>& factors) {
    if (factors.empty()) return make_unique<Constant>(1.0);
    if (factors.size() == 1) return factors[0]->simplify();
    dExp acc = make_unique<Multiply>(factors[0], factors[1]);
    for (size_t i = 2; i < factors.size(); ++i) {
        acc = make_unique<Multiply>(toShared(move(acc)), factors[i]);
    }
    return acc;
}

static bool extractCommonFactor(vector<shared_ptr<Exp>>& a,
                                vector<shared_ptr<Exp>>& b,
                                shared_ptr<Exp>& common) {
    for (size_t i = 0; i < a.size(); ++i) {
        const string aStr = a[i]->toString();
        for (size_t j = 0; j < b.size(); ++j) {
            if (aStr == b[j]->toString()) {
                common = a[i];
                a.erase(a.begin() + static_cast<long long>(i));
                b.erase(b.begin() + static_cast<long long>(j));
                return true;
            }
        }
    }
    return false;
}
static bool extractVariableFromPower(vector<shared_ptr<Exp>>& a,
                                     vector<shared_ptr<Exp>>& b,
                                     shared_ptr<Exp>& common) {
    for (size_t i = 0; i < a.size(); ++i) {
        if (!dynamic_cast<VariableX*>(a[i].get())) continue;
        for (size_t j = 0; j < b.size(); ++j) {
            auto p = dynamic_cast<Power*>(b[j].get());
            if (!p || p->hasFraction) continue;
            if (!isInt(p->exponent) || p->exponent < 1) continue;

            common = make_shared<VariableX>();
            double nextExp = p->exponent - 1;
            a.erase(a.begin() + static_cast<long long>(i));
            b.erase(b.begin() + static_cast<long long>(j));
            if (nextExp > 0) {
                if (nextExp == 1) b.insert(b.begin(), make_shared<VariableX>());
                else b.insert(b.begin(), make_shared<Power>(nextExp));
            }
            return true;
        }
    }
    return false;
}

struct Poly {
    bool ok = true;
    map<int, double> terms;
};

static Poly polyAdd(const Poly& a, const Poly& b, double sign) {
    Poly out;
    if (!a.ok || !b.ok) {
        out.ok = false;
        return out;
    }
    out.terms = a.terms;
    for (const auto& kv : b.terms) {
        out.terms[kv.first] += sign * kv.second;
    }
    return out;
}

static Poly polyMul(const Poly& a, const Poly& b) {
    Poly out;
    if (!a.ok || !b.ok) {
        out.ok = false;
        return out;
    }
    for (const auto& ka : a.terms) {
        for (const auto& kb : b.terms) {
            out.terms[ka.first + kb.first] += ka.second * kb.second;
        }
    }
    return out;
}

static Poly toPoly(const Exp* expr) {
    Poly out;
    if (auto c = dynamic_cast<const Constant*>(expr)) {
        out.terms[0] = c->value;
        return out;
    }
    if (dynamic_cast<const VariableX*>(expr)) {
        out.terms[1] = 1.0;
        return out;
    }
    if (auto p = dynamic_cast<const Power*>(expr)) {
        if (isInt(p->exponent) && p->exponent >= 0) {
            out.terms[static_cast<int>(llround(p->exponent))] = 1.0;
            return out;
        }
        out.ok = false;
        return out;
    }
    if (auto add = dynamic_cast<const AddSub*>(expr)) {
        auto l = toPoly(add->left.get());
        auto r = toPoly(add->right.get());
        return polyAdd(l, r, add->op == '+' ? 1.0 : -1.0);
    }
    if (auto mul = dynamic_cast<const Multiply*>(expr)) {
        auto l = toPoly(mul->left.get());
        auto r = toPoly(mul->right.get());
        return polyMul(l, r);
    }
    out.ok = false;
    return out;
}

static Poly toPoly(const shared_ptr<Exp>& expr) {
    return toPoly(expr.get());
}

static dExp polyToExpr(const Poly& p) {
    dExp acc;
    for (auto it = p.terms.rbegin(); it != p.terms.rend(); ++it) {
        double coeff = it->second;
        int exp = it->first;
        if (fabs(coeff) < 1e-12) continue;
        bool negative = coeff < 0.0;
        double abscoeff = fabs(coeff);

        dExp term;
        if (exp == 0) {
            term = make_unique<Constant>(abscoeff);
        } else {
            dExp base;
            if (exp == 1) base = make_unique<VariableX>();
            else base = make_unique<Power>(exp);
            if (abscoeff == 1.0) term = move(base);
            else term = make_unique<Multiply>(make_shared<Constant>(abscoeff), toShared(move(base)));
        }

        if (!acc) {
            if (negative) {
                if (exp == 0) acc = make_unique<Constant>(-abscoeff);
                else acc = make_unique<Multiply>(make_shared<Constant>(-1.0), toShared(move(term)));
            } else {
                acc = move(term);
            }
        } else {
            acc = make_unique<AddSub>(
                toShared(move(acc)),
                toShared(move(term)),
                negative ? '-' : '+'
            );
        }
    }

    if (!acc) return make_unique<Constant>(0.0);
    return acc;
}

Constant::Constant(double v) : value(v) {}
Constant::Constant(long long n, long long d) : value(static_cast<double>(n) / static_cast<double>(d)) {
    hasFraction = true;
    num = n;
    den = d;
    normaliseFraction(num, den);
}
string Constant::toString() const {
    if (hasFraction) {
        return formatFraction(num, den);
    }
    return formatNumber(value);
}
dExp Constant::derivative() const {
    return make_unique<Constant>(0);
}
dExp Constant::simplify() const {
    if (hasFraction) return make_unique<Constant>(num, den);
    return make_unique<Constant>(value);
}
double Constant::evaluate(double x) const {
    return value;
}

string VariableX::toString() const {
    return "x";
}
dExp VariableX::derivative() const {
    return make_unique<Constant>(1);
}
dExp VariableX::simplify() const {
    return make_unique<VariableX>();
}
double VariableX::evaluate(double x) const {
    return x;
}

Power::Power(double n) : exponent(n) {}
Power::Power(long long n, long long d) : exponent(static_cast<double>(n) / static_cast<double>(d)) {
    hasFraction = true;
    num = n;
    den = d;
    normaliseFraction(num, den);
}
string Power::toString() const {
    if (hasFraction) {
        if (den == 1) return "x^" + to_string(num);
        return "x^(" + formatFraction(num, den) + ")";
    }
    return "x^" + formatNumber(exponent);
}
dExp Power::derivative() const {
    if (hasFraction) {
        long long n = num;
        long long d = den;
        long long n_minus = n - d;
        return make_unique<Multiply>(
            make_unique<Constant>(n, d),
            make_unique<Power>(n_minus, d)
        )->simplify();
    }
    return make_unique<Multiply>(
        make_unique<Constant>(exponent),
        make_unique<Power>(exponent - 1)
    )->simplify();
}
dExp Power::simplify() const {
    if (hasFraction) {
        if (num == 0) return make_unique<Constant>(1);
        if (den == 1 && num == 1) return make_unique<VariableX>();
        return make_unique<Power>(num, den);
    }
    if (exponent == 0) return make_unique<Constant>(1);
    if (exponent == 1) return make_unique<VariableX>();
    return make_unique<Power>(exponent);
}
double Power::evaluate(double x) const {
    return pow(x, exponent);
}

Exponential::Exponential(double a) : coefficient(a) {}
string Exponential::toString() const {
    if (coefficient == 1) return "e^x";
    return "e^(" + formatNumber(coefficient) + "*x)";
}
dExp Exponential::derivative() const {
    return make_unique<Multiply>(
        make_unique<Constant>(coefficient),
        make_unique<Exponential>(coefficient)
    )->simplify();
}
dExp Exponential::simplify() const {
    return make_unique<Exponential>(coefficient);
}
double Exponential::evaluate(double x) const {
    return exp(coefficient * x);
}

AddSub::AddSub(shared_ptr<Exp> l, shared_ptr<Exp> r, char o) : left(l), right(r), op(o) {}
string AddSub::toString() const {
    return left->toString() + " " + op + " " + right->toString();
}
dExp AddSub::derivative() const {
    return make_unique<AddSub>(left->derivative(), right->derivative(), op)->simplify();
}
dExp AddSub::simplify() const {
    auto l = left->simplify();
    auto r = right->simplify();
    auto p = polyAdd(toPoly(l.get()), toPoly(r.get()), op == '+' ? 1.0 : -1.0);
    if (p.ok) return polyToExpr(p);

    auto lShared = toShared(move(l));
    auto rShared = toShared(move(r));
    auto lc = asConst(lShared);
    auto rc = asConst(rShared);

    if (lc && rc) {
        long long ln, ld, rn, rd;
        if (getRational(lc, ln, ld) && getRational(rc, rn, rd)) {
            long long num = (op == '+') ? (ln * rd + rn * ld) : (ln * rd - rn * ld);
            long long den = ld * rd;
            return makeRationalConst(num, den);
        }
        double v = (op == '+') ? (lc->value + rc->value) : (lc->value - rc->value);
        return make_unique<Constant>(v);
    }

    if (op == '+') {
        if (lc && lc->value == 0.0) return move(r);
        if (rc && rc->value == 0.0) return move(l);
    } else {
        if (rc && rc->value == 0.0) return move(l);
        if (lc && lc->value == 0.0) {
            if (rc) {
                long long rn, rd;
                if (getRational(rc, rn, rd)) {
                    return makeRationalConst(-rn, rd);
                }
                return make_unique<Constant>(-rc->value);
            }
            return make_unique<Multiply>(make_shared<Constant>(-1.0), rShared)->simplify();
        }
    }

    if (op == '-' && isTanTimesOnePlusTan(lShared.get()) && isSecSquaredExpr(rShared.get())) {
        return make_unique<AddSub>(
            make_unique<Tangent>(),
            make_unique<Constant>(1),
            '-'
        )->simplify();
    }
    if (op == '-' && isSecSquaredExpr(lShared.get()) && isTanTimesOnePlusTan(rShared.get())) {
        return make_unique<AddSub>(
            make_unique<Constant>(1),
            make_unique<Tangent>(),
            '-'
        )->simplify();
    }

    vector<shared_ptr<Exp>> lf;
    vector<shared_ptr<Exp>> rf;
    collectFactors(lShared, lf);
    collectFactors(rShared, rf);

    shared_ptr<Exp> common;
    if (extractCommonFactor(lf, rf, common)) {
        if (!isConstValue(common.get(), 1.0) && !isConstValue(common.get(), -1.0)) {
            auto restL = buildProduct(lf);
            auto restR = buildProduct(rf);
            auto inner = make_shared<AddSub>(restL, restR, op);
            return make_unique<Multiply>(common, inner)->simplify();
        }
    }
    if (extractVariableFromPower(lf, rf, common)) {
        auto restL = buildProduct(lf);
        auto restR = buildProduct(rf);
        auto inner = make_shared<AddSub>(restL, restR, op);
        return make_unique<Multiply>(common, inner)->simplify();
    }
    if (extractVariableFromPower(rf, lf, common)) {
        auto restL = buildProduct(rf);
        auto restR = buildProduct(lf);
        auto inner = make_shared<AddSub>(restL, restR, op);
        return make_unique<Multiply>(common, inner)->simplify();
    }

    return make_unique<AddSub>(lShared, rShared, op);
}
double AddSub::evaluate(double x) const {
    if (op == '+') return left->evaluate(x) + right->evaluate(x);
    return left->evaluate(x) - right->evaluate(x);
}

Multiply::Multiply(shared_ptr<Exp> l, shared_ptr<Exp> r) : left(l), right(r) {}
string Multiply::toString() const {
    vector<shared_ptr<Exp>> factors;
    collectFactors(left, factors);
    collectFactors(right, factors);

    if (factors.size() == 2) {
        if (dynamic_cast<Sine*>(factors[0].get()) && dynamic_cast<Sine*>(factors[1].get())) return "sin^2(x)";
        if (dynamic_cast<Cosine*>(factors[0].get()) && dynamic_cast<Cosine*>(factors[1].get())) return "cos^2(x)";
        if (dynamic_cast<Tangent*>(factors[0].get()) && dynamic_cast<Tangent*>(factors[1].get())) return "tan^2(x)";
        if (dynamic_cast<Secant*>(factors[0].get()) && dynamic_cast<Secant*>(factors[1].get())) return "sec^2(x)";
        if (dynamic_cast<Cosecant*>(factors[0].get()) && dynamic_cast<Cosecant*>(factors[1].get())) return "csc^2(x)";
        if (dynamic_cast<Cotangent*>(factors[0].get()) && dynamic_cast<Cotangent*>(factors[1].get())) return "cot^2(x)";
    }

    bool hasTrigOrExp = false;
    for (const auto& f : factors) {
        if (isTrigLikeExpr(f.get()) || isExponentialLikeExpr(f.get())) {
            hasTrigOrExp = true;
            break;
        }
    }

    vector<shared_ptr<Exp>> consts;
    vector<shared_ptr<Exp>> addsubs;
    vector<shared_ptr<Exp>> others;
    for (const auto& f : factors) {
        if (asConst(f)) {
            consts.push_back(f);
        } else if (isAddSubOrDivideExpr(f.get())) {
            addsubs.push_back(f);
        } else if (hasTrigOrExp && isPowerLikeExpr(f.get())) {
            addsubs.push_back(f);
        } else {
            others.push_back(f);
        }
    }

    dExp constProd;
    if (!consts.empty()) {
        constProd = buildProductUnique(consts);
        constProd = constProd->simplify();
        if (auto c = asConst(toShared(move(constProd)))) {
            if (c->value == 0.0) return "0";
            if (c->value == 1.0) {
                constProd.reset();
            }
        }
    }

    vector<string> parts;
    if (constProd) parts.push_back(constProd->toString());
    for (const auto& o : others) parts.push_back(o->toString());
    for (const auto& a : addsubs) parts.push_back("(" + a->toString() + ")");

    if (parts.empty()) return "1";
    string out = parts[0];
    for (size_t i = 1; i < parts.size(); ++i) {
        out += "*" + parts[i];
    }
    return out;
}
dExp Multiply::derivative() const {
    return make_unique<AddSub>(
        make_unique<Multiply>(left->derivative(), right),
        make_unique<Multiply>(left, right->derivative()),
        '+'
    )->simplify();
}
dExp Multiply::simplify() const {
    auto l = left->simplify();
    auto r = right->simplify();

    auto lShared = toShared(move(l));
    auto rShared = toShared(move(r));
    auto lc = asConst(lShared);
    auto rc = asConst(rShared);

    if (lc && rc) {
        long long ln, ld, rn, rd;
        if (getRational(lc, ln, ld) && getRational(rc, rn, rd)) {
            return makeRationalConst(ln * rn, ld * rd);
        }
        return make_unique<Constant>(lc->value * rc->value);
    }
    if (lc && lc->value == 0.0) return make_unique<Constant>(0);
    if (rc && rc->value == 0.0) return make_unique<Constant>(0);
    if (lc && lc->value == 1.0) return rShared->simplify();
    if (rc && rc->value == 1.0) return lShared->simplify();
    if (lc && lc->value == -1.0) {
        return make_unique<Multiply>(make_shared<Constant>(-1.0), rShared)->simplify();
    }
    if (rc && rc->value == -1.0) {
        return make_unique<Multiply>(make_shared<Constant>(-1.0), lShared)->simplify();
    }

    vector<shared_ptr<Exp>> lf;
    vector<shared_ptr<Exp>> rf;
    collectFactors(lShared, lf);
    collectFactors(rShared, rf);

    vector<shared_ptr<Exp>> all = lf;
    all.insert(all.end(), rf.begin(), rf.end());

    vector<shared_ptr<Exp>> consts;
    vector<shared_ptr<Exp>> nonconsts;
    for (auto& f : all) {
        if (asConst(f)) consts.push_back(f);
        else nonconsts.push_back(f);
    }

    dExp constProd;
    if (!consts.empty()) {
        constProd = buildProductUnique(consts);
        constProd = constProd->simplify();
        if (auto c = asConst(toShared(move(constProd)))) {
            if (c->value == 0.0) return make_unique<Constant>(0);
            if (c->value == 1.0) constProd.reset();
        }
    }

    vector<shared_ptr<Exp>> merged;
    if (constProd) merged.push_back(toShared(move(constProd)));
    for (auto& f : nonconsts) merged.push_back(f);

    return buildProductUnique(merged);
}
double Multiply::evaluate(double x) const {
    return left->evaluate(x) * right->evaluate(x);
}

Divide::Divide(shared_ptr<Exp> l, shared_ptr<Exp> r) : left(l), right(r) {}
string Divide::toString() const {
    return "(" + left->toString() + ")/(" + right->toString() + ")";
}
dExp Divide::derivative() const {
    return make_unique<Divide>(
        make_unique<AddSub>(
            make_unique<Multiply>(left->derivative(), right),
            make_unique<Multiply>(left, right->derivative()),
            '-'
        ),
        make_unique<Multiply>(right, right)
    )->simplify();
}
dExp Divide::simplify() const {
    auto l = left->simplify();
    auto r = right->simplify();
    auto lShared = toShared(move(l));
    auto rShared = toShared(move(r));
    auto lc = asConst(lShared);
    auto rc = asConst(rShared);

    if (lc && rc) {
        long long ln, ld, rn, rd;
        if (getRational(lc, ln, ld) && getRational(rc, rn, rd) && rn != 0) {
            return makeRationalConst(ln * rd, ld * rn);
        }
        return make_unique<Constant>(lc->value / rc->value);
    }
    if (lc && lc->value == 0.0) return make_unique<Constant>(0);
    if (rc && rc->value == 1.0) return move(lShared)->simplify();
    if (rc && rc->value == -1.0) {
        return make_unique<Multiply>(make_shared<Constant>(-1.0), lShared)->simplify();
    }

    return make_unique<Divide>(lShared, rShared);
}
double Divide::evaluate(double x) const {
    double denom = right->evaluate(x);
    if (denom == 0) return NAN;
    return left->evaluate(x) / denom;
}

dExp Constant::substitute(const shared_ptr<Exp>& replacement) const {
    if (hasFraction) return make_unique<Constant>(num, den);
    return make_unique<Constant>(value);
}
dExp VariableX::substitute(const shared_ptr<Exp>& replacement) const {
    return dExp(move(replacement->simplify()));
}
dExp Power::substitute(const shared_ptr<Exp>& replacement) const {
    if (hasFraction) {
        return make_unique<PowerComposed>(replacement, num, den)->simplify();
    }
    return make_unique<PowerComposed>(replacement, exponent)->simplify();
}
dExp Exponential::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<Exponential>(coefficient);
}
dExp AddSub::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<AddSub>(
        left->substitute(replacement),
        right->substitute(replacement),
        op
    )->simplify();
}
dExp Multiply::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<Multiply>(
        left->substitute(replacement),
        right->substitute(replacement)
    )->simplify();
}
dExp Divide::substitute(const shared_ptr<Exp>& replacement) const {
    return make_unique<Divide>(
        left->substitute(replacement),
        right->substitute(replacement)
    )->simplify();
}

#endif