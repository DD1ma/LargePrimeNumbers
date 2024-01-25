#pragma once

namespace BPSW {

template <class T>
int chooseD(T num);

template <class T>
struct LucasSequencesValues {
    int D;
    T U, V, Q;
};

template <class T>
LucasSequencesValues<T> get_method_A_values(T n);

template <class T>
LucasSequencesValues<T> get_method_A_star_values(T n);

template <class T>
LucasSequencesValues<T> LucasAStarTest(T &n, LucasSequencesValues<T> initial_UVQ, T &MOD);

template <class T>
bool BailliePSWTest(T num);

};  // namespace BPSW
