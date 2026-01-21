#ifndef EXPRESSION_HPP
#define EXPRESSION_HPP

#include <memory>
#include <string>

using namespace std;

class Exp {
    public:
        virtual ~Exp() = default;
        virtual string toString() const = 0;
        virtual unique_ptr<Exp> derivative() const = 0;
        virtual unique_ptr<Exp> simplify() const = 0;
        virtual double evaluate(double x) const = 0;
        virtual unique_ptr<Exp> substitute(const shared_ptr<Exp>& replacement) const = 0;
};

using dExp = unique_ptr<Exp>;


#endif
