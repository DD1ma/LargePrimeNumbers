#include "allchecks.h"
#include "Baillie-PSW.h"
#include <vector>

namespace BPSW {

template <class T>
int chooseD(T num) {
    int D = 5;
    int cnt = 0;
    // on average takes about 3.1477 checks
    while (MyFunctions::JacobiSymbol<T>(D, num) != -1) {
        if (D > 0) {
            D += 2;
        } else {
            D -= 2;
        }
        D *= -1;
        cnt++;
        // prevent infinite loop if n is square
        if (cnt == 5) {
            auto x = MyFunctions::GetSqrt(num);
            if (num == x * x) {
                return 0;
            }
        }
    }
    return D;
}

template <class T>
LucasSequencesValues<T> get_method_A_values(T n) {
    LucasSequencesValues<T> starting_values;
    starting_values.D = chooseD<T>(n);
    starting_values.U = 1;
    starting_values.V = 1;  // P
    starting_values.Q = (1 - D) / 4;
    return starting_values;
}

template <class T>
LucasSequencesValues<T> get_method_A_star_values(T n) {
    LucasSequencesValues<T> starting_values;
    starting_values.D = chooseD<T>(n);
    if (D == 5) {
        starting_values.U = 1;
        starting_values.V = 5;  // P
        starting_values.Q = 5;
    } else {
        starting_values.U = 1;
        starting_values.V = 1;  // P
        starting_values.Q = (1 - D) / 4;
    }
    return starting_values;
}

template <class T>
LucasSequencesValues<T> LucasAStarTest(T &n, LucasSequencesValues<T> initial_UVQ, T &MOD) {
    // we want to find U_(n + 1) % n
    if (n == 1) {
        return initial_UVQ;
    }
    // It should be better to just add MOD instead of multiplying by 1/2
    T one_half = MyFunctions::binpow<T>(2, MOD - 2, MOD);
    LucasSequencesValues<T> UVQ_now;
    if (n % 2 == 0) {
        // U_2k = U_k * V_k
        // V_2k = ((V_k)^2 + D(U_k)^2) / 2
        // Q_2k = (Q_k)^2
        LucasSequencesValues<T> UVQ_half = LucasAStarTest(n / 2, D, P, U, V, MOD);
        UVQ_now.U = UVQ_half.U * UVQ_half.V % MOD;
        UVQ_now.V = (UVQ_half.V * UVQ_half.V + D * UVQ_half.U * UVQ_half.U) % MOD * one_half % MOD;
        UVQ_now.Q = UVQ_half.Q * UVQ_half.Q % MOD;
        return UVQ_now;
    }
    // U_(2k+1) = (P * U_2k + V_2k) / 2
    // V_(2k+1) = (D * U_2k + P * V_2k) / 2
    // Q_(2k+1) = Q_2k * Q
    LucasSequencesValues<T> UVQ_prev = LucasAStarTest(n - 1, D, P, U, V, MOD);
    UVQ_now.U = (initial_UVQ.V * UVQ_prev.U + UVQ_prev.V) % MOD * one_half % MOD;
    UVQ_now.V = (D * UVQ_prev.U + initial_UVQ.V * UVQ_prev.V) % MOD * one_half % MOD;
    UVQ_now.Q = UVQ_prev.Q * initial_UVQ.Q % MOD;
    return UVQ_now;
}

template <class T>
bool BailliePSWTest(T num) {
    // 1. Trial factoring
    const std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 43};
    for (auto prime : primes) {
        if (num % prime == 0) {
            return false;
        }
    }
    // 2. MillerRabbin base 2
    if (!Checks::MillerRabbin()) {
        return false;
    }
    // 3. find D
    LucasSequencesValues<T> UVQ;
    UVQ = get_method_A_star_values(num);
    if (UVQ.D == 0) {
        return false;
    }
    // 4. Lucas Strong test
    LucasSequencesValues<T> final_UVQ = LucasAStarTest<T>(num + 1, UVQ, num);
    if (final_UVQ.U != 0) {
        return false;
    }
    return true;
}
};  // namespace BPSW
