#pragma once

namespace Checks {

bool mersenne_trial_factoring(int p);

template <class T>
inline bool check_mersenne(int p);

template <class T>
inline bool PRP_LL(int p);

template <class T>
bool MillerRabbin(T &p, int tests);

template <class T>
inline bool FermatPrimalityTest(T &p);

template <class T>
bool PollardsFactorization(int p);

};  // namespace Checks
