#include "chain_rule.cpp"
#include "polynomials_and_exponential_functions.cpp"
#include "trigonometric_functions.cpp"
#include "inverse_trigonometric_functions.cpp"
#include "implicit_differentiation.cpp"
#include "expression_utils.hpp"

#ifndef MAIN_CPP
#define MAIN_CPP

#include <iostream>
using namespace std;

int main() {
    // Find y' if sin(x + y) = (y^2) * cos(x).
    ImplicitEquation equation = ImplicitEquation(
        make_unique<SineComposed>(
            make_unique<AddSub>(
                make_shared<VariableX>(),
                make_shared<VariableY>(),
                '+'
            )
        ),
        make_unique<Multiply>(
            make_unique<Power>(make_shared<VariableY>(), 2),
            make_unique<CosineComposed>(make_shared<VariableX>())
        )
    );
    cout << "The implicit equation is: " << equation.toString() << endl;
    dExp derivative = equation.derivative();
    cout << "The derivative dy/dx is: " << derivative->toString() << endl;

    return 0;
}

#endif