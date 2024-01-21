// c++ -I boost_1_83_0 main.cpp -o main -lgmp -O2
// time ./main
#include "tests.h"
#include "sieve.h"

// TODO make every check in its own program
// TODO multithread
// TODO just casual speed up for basic functions.

int main() {
    const int Max_Precalc = 1e6;
    sieve::PrecalcPrimes(Max_Precalc);
    Tests::run_all_tests();
}
