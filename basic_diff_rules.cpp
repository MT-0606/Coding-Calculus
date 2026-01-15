#include "basic_diff_rules.hpp"

/***************************************
 * RULE 1: (c)' = 0
 ***************************************/
Constant::Constant(double v) : value(v) {
    cout << "Computing the derivative of " << value << endl;
}
string Constant::toString() const {
    return formatNumber(value);
}
dExp Constant::derivative() const {
    return make_unique<Constant>(0);
}
dExp Constant::simplify() const {
    return make_unique<Constant>(value);
}
double Constant::evaluate(double x) const {
    return value;
}

/***************************************
 * RULE 2: (x)' = 1
 ***************************************/
string Variable::toString() const {
    return "x";
}
dExp Variable::derivative() const {
    return make_unique<Constant>(1);
}
dExp Variable::simplify() const {
    return make_unique<Variable>();
}
double Variable::evaluate(double x) const {
    return x;
}

/***************************************
 * RULE 3: (x^n)' = n*x^(n-1)
 ***************************************/
Power::Power(double n) : exponent(n) {
    cout << "Computing the derivative of x^" << exponent << endl;
}
string Power::toString() const {
    return "x^" + formatNumber(exponent);
}
dExp Power::derivative() const {
    // (x^n)' = n*x^(n-1)
    return make_unique<Multiply>(
        make_unique<Constant>(exponent),
        make_unique<Power>(exponent - 1)
    );
}
dExp Power::simplify() const {
    if (exponent == 0) return make_unique<Constant>(1); // x^0 = 1
    if (exponent == 1) return make_unique<Variable>();  // x^1 = x
    return make_unique<Power>(exponent);
}
double Power::evaluate(double x) const {
    return std::pow(x, exponent);
}

/***************************************
 * RULE 4: (e^(a*x))' = a*e^(a*x)
 ***************************************/
Exponential::Exponential(double a) : coefficient(a) {
    cout << "Computing the derivative of e^(" << coefficient << "*x)" << endl;
}
string Exponential::toString() const {
    if (coefficient == 1) return "e^x";
    return "e^(" + formatNumber(coefficient) + "*x)";
}
dExp Exponential::derivative() const {
    // (e^(a*x))' = a*e^(a*x)
    return make_unique<Multiply>(
        make_unique<Constant>(coefficient),
        make_unique<Exponential>(coefficient)
    );
}
dExp Exponential::simplify() const {
    return make_unique<Exponential>(coefficient);
}

double Exponential::evaluate(double x) const {
    return std::exp(coefficient * x);
}
/***************************************
 * RULE 5: (f±g)' = f'±g'
 ***************************************/
AddSub::AddSub(shared_ptr<Exp> l, shared_ptr<Exp> r, char o) : left(l), right(r), op(o) {
    cout << "Computing the derivative of (" << left->toString() << " " << op << " " << right->toString() << ")" << endl;
}
string AddSub::toString() const {
    return "(" + left->toString() + " " + op + " " + right->toString() + ")";
}
dExp AddSub::derivative() const {
    return make_unique<AddSub>(left->derivative(), right->derivative(), op);
}
dExp AddSub::simplify() const {
    auto l = left->simplify();
    auto r = right->simplify();
    return make_unique<AddSub>(l, r, op);
}
double AddSub::evaluate(double x) const {
    if (op == '+') return left->evaluate(x) + right->evaluate(x);
    return left->evaluate(x) - right->evaluate(x);  // op == '-'
}

/***************************************
 * RULE 6: (f*g)' = f'*g + f*g'
 ***************************************/
Multiply::Multiply(shared_ptr<Exp> l, shared_ptr<Exp> r) : left(l), right(r) {
    cout << "Computing the derivative of (" << left->toString() << " * " << right->toString() << ")" << endl;
}
string Multiply::toString() const {
    return "(" + left->toString() + " * " + right->toString() + ")";
}
dExp Multiply::derivative() const {
    return make_unique<AddSub>(
        make_unique<Multiply>(left->derivative(), right),
        make_unique<Multiply>(left, right->derivative()),
        '+'
    );
}
dExp Multiply::simplify() const {
    auto l = left->simplify();
    auto r = right->simplify();
    return make_unique<Multiply>(l, r);
}
double Multiply::evaluate(double x) const {
    return left->evaluate(x) * right->evaluate(x);
}

/***************************************
 * RULE 7: (f/g)' = (f'*g - f*g')/(g^2)
 ***************************************/
Divide::Divide(shared_ptr<Exp> l, shared_ptr<Exp> r) : left(l), right(r) {
    cout << "Computing the derivative of (" << left->toString() << " / " << right->toString() << ")" << endl;
}
string Divide::toString() const {
    return "(" + left->toString() + " / " + right->toString() + ")";
}
dExp Divide::derivative() const {
    return make_unique<Divide>(
        make_unique<AddSub>(
            make_unique<Multiply>(left->derivative(), right),
            make_unique<Multiply>(left, right->derivative()),
            '-' // numerator: f'*g - f*g'
        ),
        make_unique<Power>(2)  // denominator: g^2
    );
}
dExp Divide::simplify() const {
    auto l = left->simplify();
    auto r = right->simplify();
    return make_unique<Divide>(l, r);
}
double Divide::evaluate(double x) const {
    double denom = right->evaluate(x);
    if (denom == 0) return NAN;
    return left->evaluate(x) / denom;
}