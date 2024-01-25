#include "utils.h"
#include <cassert>
#include <cmath>

namespace MyFunctions {

template <class T>
inline T binpow(T &a, T &b, T &MOD) {
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
bool GetSqrt(T &number) {
    assert(number >= 0);  // negative number cannot be a square
    // Heron`s method (Newton`s method for x^2 - number)
    // x_{n + 1} = 1/2(x_n + number / x_n)
    return std::sqrt(number);
}

// TODO faster remainder by modulo 2^p - 1
// And check its speed with %

template <class T>
inline T Gcd(T a, T b) {
    while (b > 0) {
        T q = a / b;
        T r = a - b * q;
        a = b;
        b = r;
    }
    return a;
}

template <class T>
T get_high_order_digit(T n, T &base, int &len) {
    while (n >= base) {
        n /= base;
        len++;
    }
    return n;
}

// https://www.imsc.res.in/~kapil/crypto/notes/node11.html
template <class T>
T LehmerGCD(T x, T y, T &base) {
    if (x < y) {
        std::swap(x, y);
    }
    T temp;
    while (y >= base) {
        int xlen = 0, ylen = 0;
        T x1 = get_high_order_digit(x, base, xlen);
        T y1 = get_high_order_digit(y, base, ylen);
        T a = 1, b = 0, c = 0, d = 1;
        if (xlen == ylen) {
            while (y1 + c != 0 && y1 + d != 0) {
                T q = (x + a) / (y + c);
                T q1 = (x + b) / (y + d);
                if (q != q1) {
                    break;
                }
                temp = a; a = c; c = temp - q * c;
                temp = b; b = d; d = b - q * d;
                temp = x1; x1 = y1; y1 = x1 - q * y1;
            }
        }
        if (b == 0) {
            temp = x; x = y; y = temp % y;
        } else {
            temp = x; x = a * x + b * y; y = c * temp + d * y;
        }
    }
    return Gcd(x, y);
}

};  // namespace MyFunctions
