#pragma once

namespace MyFunctions {

template <class T>
inline T binpow(T a, T b, T &MOD);

template <class T>
int JacobiSymbol(T a, T n);

template <class T>
bool GetSqrt(T &number);

template <class T>
inline T Gcd(T a, T b);

template <class T>
T LehmerGCD(T a, T b, T &base);

};  // namespace MyFunctions
