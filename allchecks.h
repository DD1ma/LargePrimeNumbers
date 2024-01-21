#pragma once

#include <boost/multiprecision/gmp.hpp>

typename boost::multiprecision::mpz_int mpz_int

bool mersenne_trial_factoring(int p);

inline bool check_mersenne(int p);

inline bool PRP_LL(int p);

bool MillerRabbin(mpz_int p, int tests);

inline bool FermatPrimalityTest(mpz_int p);

bool PollardsFactorization(int p);
