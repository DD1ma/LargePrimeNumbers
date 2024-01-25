#include "allchecks.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include "tests.h"
#include "sieve.h"

// TODO !! class for time measuring

namespace Tests {

template <class T>
void test_target_Mersenne(int p) {
    const auto start{std::chrono::steady_clock::now()};
    if (Checks::check_mersenne<T>(p)) {
        const auto end{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_seconds{end - start};
        std::cout << p << " is Mersenne prime" << '\n';
        std::cout << elapsed_seconds.count() << '\n';
    }
}

template <class T>
void test_section_Mersenne() {
    const int Max_Mersenne = 1e5;
    std::cout << std::fixed << std::setprecision(10);
    for (int i = 3; i <= Max_Mersenne; ++i) {
        if (sieve::is_prime[i]) {
            test_target_Mersenne<T>(i);
        }
    }
    // std::cout << "Mersenne end\n";
}

void test_section_Miller_Rabin() {
    const int MaxMiller = 1e6;
    for (int i = 5; i <= MaxMiller; i += 2) {
        if (Checks::MillerRabbin<int>(i, 10) != sieve::is_prime[i]) {
            std::cout << "We Have Some Problems with p = " << i << '\n';
        }
    }
    // std::cout << "Miller-Rabin end\n";
}

// TODO Need to add assert tests

template <class T>
void run_all_tests() {
    // test_target_Mersenne();
    // test_section_Mersenne();
    // test_section_Miller_Rabin();
    test_target_Mersenne<T>(132049);  // true
    test_target_Mersenne<T>(44497);   // true
    test_target_Mersenne<T>(100003);  // false
}

};  // namespace Tests
