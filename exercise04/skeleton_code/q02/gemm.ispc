// File       : gemm.ispc
// Description: ISPC GEMM body
// Copyright 2020 ETH Zurich. All Rights Reserved.
#include "common.h"

/**
 * @brief General matrix-matrix multiplication kernel (GEMM). Computes C = AB.
 * ISPC implementation used for SSE2 and AVX2 extended ISA
 *
 * @param A Matrix dimension p x r
 * @param B Matrix dimension r x q
 * @param C Matrix dimension p x q
 * @param p Dimensional parameter
 * @param r Dimensional parameter
 * @param q Dimensional parameter
 */
/* TODO: Missing piece */
#ifdef _ISPC_SSE2_
void gemm_sse2(
#else
void gemm_avx2(
#endif /* _ISPC_SSE2_ */
        const /* TODO: Missing piece */ A,
        const /* TODO: Missing piece */ B,
        /* TODO: Missing piece */ C,
        const /* TODO: Missing piece */ int p,
        const /* TODO: Missing piece */ int r,
        const /* TODO: Missing piece */ int q)
{
    ///////////////////////////////////////////////////////////////////////////
    // TODO: Write your ISPC implementation of a matrix-matrix multiplication
    // here.  Try to use a general type for floating point numbers, such that
    // you can easily compile this kernel for float (single precision, 32bit)
    // and double (double precision, 64bit).  Check out the 'common.h' header
    // for a hint.  Note: A working code can be implemented with ~15 lines of
    // ISPC code.
    //
    // We use two different function names for the same function, such that we
    // can compile for two different targets (SSE2 and AVX2) and link both of
    // them to our application code in 'gemm.cpp'.
    //
    // Your ISPC implementation should correspond to an optimized serial 
    // algorithm (for example, one the uses blocked matrix multiplication).
    // Otherwise, the measured speedup might be unrealistically large.
    ///////////////////////////////////////////////////////////////////////////
}
