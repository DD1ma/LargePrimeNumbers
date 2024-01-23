#include "utils.h"
#include <cassert>
#include <cmath>

namespace MyFunctions {

template <class T>
inline T binpow(T a, T b, T &MOD) {
    assert(b >= 0 && "Power must be at least zero");
    assert(MOD >= 1 && "Module must be at least one");
    T ans = 1;
    while (b) {
        if (b & 1) {
            ans = ans * a % MOD;
        }
        // overflows for int and long long
        a = a * a % MOD;
        b >>= 1;
    }
    return ans;
}

template <class T>
int JacobiSymbol(T a, T n) {
    assert(a > 0 && n % 2 == 1);  // invalid value encountered in Jacobi
    a %= n;
    int t = 1;
    int r;
    while (a != 0) {
        while (a % 2 == 0) {
            a /= 2;
            r = n % 8;
            if (r == 3 || r == 5) {
                t *= -1;
            }
        }
        std::swap(a, n);
        if (a % 4 == 3 && n % 4 == 3) {
            t *= -1;
        }
        a %= n;
    }
    if (n == 1) {
        return t;
    } else {
        return 0;
    }
}

template <class T>
bool GetSqrt(T number) {
    assert(number >= 0);  // negative number cannot be a square
    // Heron`s method (Newton`s method for x^2 - number)
    // x_{n + 1} = 1/2(x_n + number / x_n)
    return std::sqrt(number);
}

// TODO faster remainder by modulo 2^p - 1
// And check its speed with %

// TODO fast GCD using Lehmer's algorithm

};  // namespace MyFunctions
