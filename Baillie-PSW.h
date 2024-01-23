#pragma once

namespace BPSW {

template <class T>
int chooseD(T num);

template <class T>
struct LucasSequencesValues {
    T U, V, Q;
};

template <class T>
LucasSequencesValues<T> getAstarUVQ(int D);

template <class T>
LucasSequencesValues<T> LucasAStarTest(T &n, int D, LucasSequencesValues<T> initial_UVQ, T &MOD);

template <class T>
bool BailliePSWTest(T num);

};  // namespace BPSW
