#include "basic_diff_rules.cpp"

int main() {
    /* Example 1: if f(x) = x^6, what is f'(x)? */
    shared_ptr<Exp> f1 = make_shared<Power>(6);
    cout << "f(x) = " << f1->toString() << endl;
    string f1_prime = f1->derivative()->toString();
    cout << "f'(x) = " << f1_prime << endl << endl;

    /* Example 2: if f(x) = 3*x^4, what is f'(x)? */
    shared_ptr<Exp> f2 = make_shared<Multiply>(
        make_shared<Constant>(3),
        make_shared<Power>(4)
    );
    cout << "f(x) = " << f2->toString() << endl;
    string f2_prime = f2->derivative()->toString();
    cout << "f'(x) = " << f2_prime << endl;


    return 0;
}   