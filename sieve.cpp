#include "sieve.h"

namespace sieve {

void PrecalcPrimes(int n) {
    is_prime.assign(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; static_cast<long long>(i) * i <= n; i++) {
        if (is_prime[i]) {
            for (int j = i * i; j <= n; j += i) {
                is_prime[j] = false;
            }
        }
    }
}

// TODO Fast Eratosthenes sieve

};  // namespace sieve
