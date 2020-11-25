// File       : vectorized_reduction.cpp
// Description: Reduction with SSE/SSE2 intrinsics
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include <cassert>
#include <chrono>
#include <random>
#include <string>
#include <iostream>
#include <cstdlib>     // posix_memalign
#include <emmintrin.h> // SSE/SSE2 intrinsics header
#include <omp.h>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Helpers, YOU DO NOT NEED TO WORK ON THESE
///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Initialize array to standard normal data
 *
 * @tparam T Real type parameter
 * @param ary Array to be initialized
 * @param N Dimension of array (length)
 * @param seed Used to seed the random number generator
 */
template <typename T>
void initialize(T * const ary, const size_t N, const size_t seed=0)
{
    default_random_engine gen(seed);
    normal_distribution<T> dist(0.0, 1.0);
    for (size_t i = 0; i < N; ++i)
        ary[i] = dist(gen);
}

/**
 * @brief Reference serial reduction (addition) of array elements
 *
 * @tparam T Real type parameter
 * @param ary Input array
 * @param N Number of elements in input array
 *
 * @return Returns reduced value (scalar)
 */
template <typename T>
static inline T gold_red(const T* const ary, const size_t N)
{
    T sum = 0.0;
    for (size_t i = 0; i < N; ++i)
        sum += ary[i];
    return sum;
}

/**
 * @brief Benchmark a test kernel versus a baseline kernel
 *
 * @tparam T Real type parameter
 * @param N Number of elements to perform reduction on (length of array)
 * @param func Function pointer to test kernel
 * @param test_name String to describe additional output
 */
template <typename T>
void benchmark_serial(const size_t N, T(*func)(const T* const, const size_t), const string test_name)
{
    T * ary;
    posix_memalign((void**)&ary, 16, N*sizeof(T));
    typedef chrono::steady_clock Clock;

    // initialize data
    initialize(ary, N);

    // reference
    T res_gold = 0.0;
    gold_red(ary, N); // warm-up
    auto t1 = Clock::now();
    // collect 10 samples.  Note: we are primarily interested in the
    // performance here, not the result
    for (int i = 0; i < 10; ++i)
        res_gold += gold_red(ary, N);
    auto t2 = Clock::now();
    const double t_gold = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();

    // test
    T res = 0.0;
    (*func)(ary, N); // warm-up
    auto tt1 = Clock::now();
    // again 10 samples
    for (int i = 0; i < 10; ++i)
        res += (*func)(ary, N);
    auto tt2 = Clock::now();
    const double t = chrono::duration_cast<chrono::nanoseconds>(tt2 - tt1).count();

    // Report
    cout << test_name << ":" << endl;
    cout << "  Data type size:     " << sizeof(T) << " byte" << endl;
    cout << "  Number of elements: " << N << endl;
    cout << "  Absolute error:     " << abs(res-res_gold)  << endl;
    cout << "  Relative error:     " << abs((res-res_gold)/res_gold)  << endl;
    cout << "  Speedup:            " << t_gold/t << endl;

    // clean up
    free(ary);
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Partially complete. TODO: YOUR TASK IS TO WRITE THE MISSING CODE. Note: A
// working code can be achieved with ~35 more lines of code for all TODO's
// below.
///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Vectorized reduction kernel using SSE intrinsics (4-way SIMD)
 *
 * @param ary Input array
 * @param N Number of elements in input array
 *
 * @return Returns reduced value (scalar)
 */
static inline float sse_red_4(const float* const ary, const size_t N)
{
    ///////////////////////////////////////////////////////////////////////////
    // TODO: Write your vectorized reduction kernel here using intrinsics from
    // the xmmintrin.h header above.  The Intel intrinsics guide is a helpful
    // reference: https://software.intel.com/sites/landingpage/IntrinsicsGuide
    ///////////////////////////////////////////////////////////////////////////
    return 0.0; // the function returns the summation of all elemnts in ary
}

/**
 * @brief Vectorized reduction kernel using SSE intrinsics (2-way SIMD)
 *
 * @param ary Input array
 * @param N Number of elements in input array
 *
 * @return Returns reduced value (scalar)
 */
static inline double sse_red_2(const double* const ary, const size_t N)
{
    ///////////////////////////////////////////////////////////////////////////
    // TODO: Write your vectorized reduction kernel for doubles (2-way SIMD)
    // here.  Note this code is very similar to what you do in sse_red_4 for
    // the 4-way SIMD case (floats)
    ///////////////////////////////////////////////////////////////////////////
    return 0.0; // the function returns the summation of all elemnts in ary
}

/**
 * @brief Benchmark a test kernel versus a baseline kernel using OpenMP to
 * exploit thread level parallelism (TLP).  See benchmark_serial above for an
 * example of a serial implementation.
 *
 * @tparam T Real type parameter
 * @param N Number of elements to perform reduction on (length of array) per
 * thread
 * @param func Function pointer to test kernel
 * @param nthreads Number of threads to run the benchmark with
 * @param test_name String to describe additional output
 */
template <typename T>
void benchmark_omp(const size_t N, T(*func)(const T* const, const size_t), const size_t nthreads, const string test_name)
{
    T * ary;
    // total length of test array is nthreads*N
    posix_memalign((void**)&ary, 16, nthreads*N*sizeof(T));
    typedef chrono::steady_clock Clock;

    // initialize data
    ///////////////////////////////////////////////////////////////////////////
    // TODO: Initialize the array 'ary' using the 'initialize' function defined
    // above.
    ///////////////////////////////////////////////////////////////////////////

    // reference (sequential)
    T res_gold = gold_red(ary, nthreads*N); // warm-up
    auto t1 = Clock::now();
    // collect 10 samples.  Note: we are primarily interested in the
    // performance here, not the result
    for (int i = 0; i < 10; ++i)
        res_gold += gold_red(ary, nthreads*N);
    auto t2 = Clock::now();
    const double t_gold = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count();

    // test
    T gres = 0.0;    // result
    double gt = 0.0; // time
    ///////////////////////////////////////////////////////////////////////////
    // TODO: Write the benchmarking code for the test kernel 'func' here.  See
    // the serial implementation 'benchmark_serial' above to get an idea.
    ///////////////////////////////////////////////////////////////////////////

    // Report
    cout << test_name << ": got " << nthreads << " threads" << endl;
    cout << "  Data type size:     " << sizeof(T) << " byte" << endl;
    cout << "  Number of elements: " << nthreads*N << endl;
    cout << "  Absolute error:     " << abs(gres-res_gold)  << endl;
    cout << "  Relative error:     " << abs((gres-res_gold)/res_gold)  << endl;
    cout << "  Speedup:            " << t_gold/gt << endl;

    // clean up
    free(ary);
}

int main(void)
{
    // Two different work size we want to test our reduction on.
    constexpr size_t N0 = (1<<15);
    constexpr size_t N1 = (1<<20);

    // Get the number of threads
    int nthreads;
#pragma omp parallel
#pragma omp master
        nthreads = omp_get_num_threads();

    // Perform the benchmarks we have written above:
    // Serial:
    //   Test the vectorized reduction kernels (2-way and 4-way) using a single
    //   core
    // OMP:
    //   Test the vectorized reduction kernels (2-way and 4-way) using OpenMP
    //   to exploit TLP.
    cout << "###############################################################################" << endl;
    cout << "TESTING SIZE N = " << N0 << endl;
    // run serial tests
    benchmark_serial<float >(N0, sse_red_4, "4-way SSE (serial)");
    benchmark_serial<double>(N0, sse_red_2, "2-way SSE (serial)");
    // run concurrent tests
    benchmark_omp<float >(N0, sse_red_4, nthreads, "4-way SSE (concurrent)");
    benchmark_omp<double>(N0, sse_red_2, nthreads, "2-way SSE (concurrent)");
    cout << "###############################################################################" << endl;

    cout << "###############################################################################" << endl;
    cout << "TESTING SIZE N = " << N1 << endl;
    // run serial tests
    benchmark_serial<float >(N1, sse_red_4, "4-way SSE (serial)");
    benchmark_serial<double>(N1, sse_red_2, "2-way SSE (serial)");
    // run concurrent tests
    benchmark_omp<float >(N1, sse_red_4, nthreads, "4-way SSE (concurrent)");
    benchmark_omp<double>(N1, sse_red_2, nthreads, "2-way SSE (concurrent)");
    cout << "###############################################################################" << endl;

    return 0;
}
