// File       : gemm.cpp
// Description: General Matrix-Matrix multiplication with ISPC
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include <cassert>
#include <iostream>
#include <chrono>
#include <string>
#include <random>
#include <cassert>
#include <cstring>
#include <cmath>

#include "common.h"
#ifdef _USE_ISPC_
#include "gemm_sse2.h" // kernel signature for the ISPC SSE2 code
#include "gemm_avx2.h" // kernel signature for the ISPC AVX2 code
using namespace ispc;
#endif /* _USE_ISPC_ */
using namespace std;

/**
 * @brief General matrix-matrix multiplication kernel (GEMM). Computes C = AB.
 *
 * @tparam T Real type parameter
 * @param A Matrix dimension p x r
 * @param B Matrix dimension r x q
 * @param C Matrix dimension p x q
 * @param p Dimensional parameter
 * @param r Dimensional parameter
 * @param q Dimensional parameter
 */
template <typename T>
static void gemm_serial(const T* const A, const T* const B, T* const C,
        const int p, const int r, const int q)
{
    ///////////////////////////////////////////////////////////////////////////
    // TODO: Write your serial implementation of a matrix-matrix multiplication
    // here. The generic type T which allows to use this function for different 
    // types T (doubles, floats etc.).
    //
    // Note that a naive serial implementation might result in unrealistic
    // speedups, when comparing with an optimized vectorized ISPC code.
    // For this reason, you are encouraged to implement blocked matrix 
    // multiplication (see also slides from tutorial 1), to optimize cache 
    // usage. Using an optimized serial algorithm is something you should do
    // anyway, before deciding to improve performance further, with 
    // vectorization.
    ///////////////////////////////////////////////////////////////////////////
}

/**
 * @brief Compute the Frobenius norm between the difference of two input
 * matrices
 *
 * @tparam T Real type parameter
 * @param M Test matrix dimension p x q
 * @param truth Reference matrix dimension p x q
 * @param p Dimensional parameter
 * @param q Dimensional parameter
 *
 * @return
 */
template <typename T>
static T validate(const T* const M, const T* const truth,
        const int p, const int q)
{
    T res = 0.0;
    for (int i = 0; i < p; ++i)
        for (int j = 0; j < q; ++j)
        {
            const T diff = M[i*q + j] - truth[i*q + j];
            res += diff * diff;
        }
    return std::sqrt(res/(p*q));
}

/**
 * @brief Initialize matrices A and B to uniform random values
 *
 * @tparam T Real type parameter
 * @param A Matrix dimension p x r
 * @param B Matrix dimension r x q
 * @param p Dimensional parameter
 * @param r Dimensional parameter
 * @param q Dimensional parameter
 * @param seed Used for random number generator
 */
template <typename T>
static void initialize(T* const A, T* const B,
        const int p, const int r, const int q, const int seed=101)
{
    default_random_engine gen(seed);
    uniform_real_distribution<T> dist(0.0, 1.0);
    for (int i = 0; i < p; ++i)
        for (int k = 0; k < r; ++k)
            A[i*r + k] = dist(gen);

    for (int k = 0; k < r; ++k)
        for (int j = 0; j < q; ++j)
            B[k*q + j] = dist(gen);
}

/**
 * @brief Benchmark a test kernel versus a baseline kernel
 *
 * @tparam T Real type parameter
 * @param p Dimensional parameter
 * @param r Dimensional parameter
 * @param q Dimensional parameter
 * @param func Function pointer to test kernel
 * @param test_name String to describe additional output
 */
template <typename T>
void benchmark(const int p, const int r, const int q,
        void (*func)(const T* const, const T* const, T* const, const int, const int, const int),
        const string test_name)
{
    T *A, *B, *C, *truth;
    posix_memalign((void**)&A, 32, p*r*sizeof(T));
    posix_memalign((void**)&B, 32, r*q*sizeof(T));
    posix_memalign((void**)&C, 32, p*q*sizeof(T));
    posix_memalign((void**)&truth, 32, p*q*sizeof(T));
    typedef chrono::steady_clock Clock;

    // initialize random matrices
    initialize(A, B, p, r, q);
    memset(C, 0, p*q*sizeof(T));
    memset(truth, 0, p*q*sizeof(T));

    // compute truth
    auto t1 = Clock::now();
    gemm_serial(A, B, truth, p, r, q);
    auto t2 = Clock::now();
    const double t_gold = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();
    const T norm_truth = validate(truth, C, p, q);

    // test kernel
    auto tt1 = Clock::now();
    (*func)(A, B, C, p, r, q);
    auto tt2 = Clock::now();
    const double t = chrono::duration_cast<chrono::nanoseconds>(tt2 - tt1).count();

    // validate
    const T diff = validate(C, truth, p, q);

    cout << test_name << ":" << endl;
    cout << "  Data type size:     " << sizeof(T) << " byte" << endl;
    cout << "  Number of elements: A=" << p*r << "; B=" << r*q << "; C=" << p*q << endl;
    cout << "  Norm of truth:      " << norm_truth << endl;
    cout << "  Error:              " << diff << endl;
    cout << "  Speedup:            " << t_gold/t << endl;

    // clean up
    free(A);
    free(B);
    free(C);
    free(truth);
}


int main(void)
{
    // problem size:
    // Matrix $A \in\mathbb{R}^{ p \times r }$
    // Matrix $B \in\mathbb{R}^{ r \times q }$
    // Matrix $C \in\mathbb{R}^{ p \times q }$
    //
    // We solve C = A * B
    constexpr int p = 512;
    constexpr int r = 1024;
    constexpr int q = 1024;

    // benchmark serial case (proof of concept, expected speedup is 1.0)
    benchmark<Real>(p, r, q, gemm_serial, "GEMM serial");

#ifdef _USE_ISPC_
    // benchmark the ISPC vectorized kernels
    benchmark<Real>(p, r, q, gemm_sse2, "GEMM ISPC SSE2");
    benchmark<Real>(p, r, q, gemm_avx2, "GEMM ISPC AVX2");
#endif /* _USE_ISPC_ */

    return 0;
}
