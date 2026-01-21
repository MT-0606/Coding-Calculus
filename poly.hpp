#include "polynomials_and_exponential_functions.hpp"
#include "expression_utils.hpp"

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

#include <map>
#include <utility>
#include <vector>
using namespace std;

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
    if (dynamic_cast<const Variable*>(expr)) {
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
            if (exp == 1) base = make_unique<Variable>();
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
