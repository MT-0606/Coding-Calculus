#ifndef EXPRESSION_UTILS_HPP
#define EXPRESSION_UTILS_HPP

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

inline string formatNumber(double v) {
    if (isnan(v)) return "nan";
    if (isinf(v)) return (v > 0) ? "inf" : "-inf";

    double intpart;
    if (modf(v, &intpart) == 0.0) {
        return to_string(static_cast<long long>(intpart));
    }

    ostringstream ss;
    ss << fixed << setprecision(8) << v;
    string s = ss.str();
    if (s.find('.') != string::npos) {
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
    }
    return s;
}

inline long long gcdll(long long a, long long b) {
    while (b != 0) {
        long long t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

inline void normaliseFraction(long long& n, long long& d) {
    if (d < 0) {
        n = -n;
        d = -d;
    }
    long long g = gcdll(n, d);
    if (g == 0) return;
    n /= g;
    d /= g;
}

inline string formatFraction(long long n, long long d) {
    normaliseFraction(n, d);
    if (d == 1) return to_string(n);
    return to_string(n) + "/" + to_string(d);
}

inline bool isIntegerDouble(double v) {
    return fabs(v - round(v)) < 1e-9;
}

inline bool isInt(double v) {
    return fabs(v - round(v)) < 1e-9;
}

#endif
