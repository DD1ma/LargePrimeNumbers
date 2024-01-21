#include "utils.h"
#include <cassert>

namespace MyFunctions {

template <class T>
inline T binpow(T a, T b, T &MOD) {
    std::assert(b >= 0 && "Power must be at least zero");
    std::assert(MOD >= 1 && "Module must be at least one");
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

// TODO faster remainder by modulo 2^p - 1
// And check its speed with %

// TODO fast GCD using Lehmer's algorithm

};  // namespace MyFunctions
