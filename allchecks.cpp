#include "allchecks.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/integer/common_factor_ct.hpp>
#include <sieve.h>
#include <vector>
#include "utils.h"

// TODO just casual speed up for basic functions.

namespace Checks {

bool mersenne_trial_factoring(int p) {
    // factoring_cost < chance_of_finding_factor * primality_test_cost
    // chance is about 1/X for finding factor from 2^X to 2^(X + 1)
    // TODO make CheckBound dependent on p
    const long long CheckBound = (1ll << 50);
    int array_size = ((CheckBound + 2 * p - 1) / (2 * p) + 63) / 64;
    const unsigned long long base_value = 0x9999999999999999;
    std::vector<unsigned long long> tocheck(array_size, base_value);
    const int MaxPrimeToCheck = 10000;
    for (int q = 3; q < MaxPrimeToCheck; q += 2) {
        if (sieve::is_prime[q]) {
            // 2p(w * 64 + bit) + 1 % q == 0 <=> w * 64 + bit = -1 / 2p
            // TODO make this bit sieve work as bit sieve
            int starting_index = q - MyFunctions::binpow<int>(2 * p, q - 2, q);
            while (starting_index < array_size * 64) {
                tocheck[starting_index / 64] &= (~(1ull << (starting_index & 63)));
                starting_index += 2 * p;
            }
        }
    }
    int reverse_bits_p = 0;
    int copyP = p;
    while (copyP) {
        reverse_bits_p = (reverse_bits_p << 1) | (copyP & 1);
        copyP >>= 1;
    }
    for (int q = 0; q < array_size; q++) {
        for (int w = 0; w < 64; w++) {
            if ((tocheck[q] >> w) & 1) {
                unsigned long long MOD = 2 * p * (64 * q + w) + 1;
                unsigned long long value = 1;
                int copyReverse = reverse_bits_p;
                while (copyReverse) {
                    // тут точно нет переполнения (есть)
                    value *= value;
                    if (copyReverse & 1) {
                        value <<= 1;
                    }
                    value %= MOD;
                    copyReverse >>= 1;
                }
                if (value == 1) {
                    return true;
                }
            }
        }
    }
    return false;
}

inline bool check_mersenne(int p) {
    mpz_int start = 4;
    mpz_int MOD = (mpz_int(1) << p) - 1;
    for (int i = 0; i < p - 2; ++i) {
        start = start * start;
        start += MOD - 2;
        while (start > MOD) {
            start = (start & MOD) + (start >> p);
        }
        // if (start > MOD) {
        //     start = (start & MOD) + (start >> p);
        //     if (start > MOD) {
        //         start = (start & MOD) + (start >> p);
        //     }
        // }
        if (start == MOD) {
            start -= MOD;
        }
    }
    if (start == 0) {
        return true;
    } else {
        return false;
    }
}

// PRP (Probable Prime)
inline bool PRP_LL(int p) {
    mpz_int value = 3;
    mpz_int MOD = (mpz_int(1) << p) - 1;
    for (int q = 0; q < p; q++) {
        value = value * value;
        while (value > MOD) {
            value = (value & MOD) + (value >> p);
        }
    }
    if (value == 3) {
        return true;
    } else {
        return false;
    }
}

bool MillerRabbin(mpz_int p, int tests) {
    // TODO write my own random number generator
    boost::random::uniform_int_distribution<mpz_int> dist(2, p - 2);
    int s = 0;
    mpz_int d = p - 1;
    while (!(d & 1)) {
        d >>= 1;
        ++s;
    }
    boost::random::mt19937 gen(228);
    for (int i = 0; i < tests; i++) {
        mpz_int a = dist(gen);
        mpz_int x = MyFunctions::binpow<mpz_int>(a, d, p);
        if (x == 1 || x == p - 1) {
            continue;
        }
        for (int j = 0; j < s; j++) {
            x = x * x % p;
            if (x == 1) {
                return false;
            }
            if (x == p - 1) {
                break;
            }
        }
        if (x != p - 1) {
            return false;
        }
    }
    return true;
}

inline bool FermatPrimalityTest(mpz_int p) {
    mpz_int value = MyFunctions::binpow<mpz_int>(2, p - 1, p);
    if (value != 1) {
        return false;
    } else {
        return true;
    }
}

// TODO Baillie-PSW primality test and its variations

bool PollardsFactorization(int p) {
    const int B1 = 10000;  // maybe 1e5 - 1e6 is better constant
    // const int B2 = 1000000;
    // I still dont know is it worth to use FFT in second part of this algorithm
    // TODO second part of Pollards p-1 factorization
    mpz_int MOD = (mpz_int(1) << p) - 1;
    mpz_int x = MyFunctions::binpow<mpz_int>(9, p, MOD);
    for (int q = 2; q < B1; q++) {
        if (sieve::is_prime[q]) {
            int pk = p;
            while (pk * p < B1) {
                pk *= p;
            }
            MyFunctions::binpow<mpz_int>(x, pk, MOD);
        }
    }
    mpz_int gcd_val = boost::integer::gcd(x - 1, MOD);
    if (gcd_val != 1) {
        return true;
    }
    return false;
}

// TODO AKS and Elliptic curves

};  // namespace Checks
