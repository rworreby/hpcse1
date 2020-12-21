#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>

const double one = 1;
extern "C" void dgemm_(const char* transa, const char* transb, const int* m,
                       const int* n, const int* k, const double* alpha,
                       const double* a, const int* lda, const double* b,
                       const int* ldb, const double* beta, double* c,
                       const int* ldc);

#ifdef HYBRID

void ompCannon(double* A, double* B, double* C, int n, int t)
{
    // TODO: implementation of the local-parallel Cannon's algorithm
    // using openMP.
}

#endif

int main(int argc, char* argv[])
{
    int myRank, rankCount;

#ifdef HYBRID
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided != MPI_THREAD_FUNNELED) {
        printf("[Error] Error initializing threaded MPI.\n");
        exit(-1);
    }
#else
    MPI_Init(&argc, &argv);
#endif

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);

    // Checking whether the matrix rankCount was passed as argument
    if (argc != 2) {
        if (myRank == 0)
            printf("[Error] Must provide matrix rankCount as parameter. E.g, "
                   "./cannon 512\n");
        MPI_Finalize();
        return 1;
    }

    // Checking that the number of MPI ranks is a square number (required by
    // Cannon's)
    int p = sqrt(rankCount);
    if (rankCount != p * p) {
        if (myRank == 0)
            printf(
                "[Error] Number of MPI Ranks must be square of an integer.\n");
        MPI_Finalize();
        return 1;
    }

    // If using MPI+OpenMP, verify that the number of threads is correct
#ifdef HYBRID
    // Getting number of threads
    int threadCount = omp_get_max_threads();

    // Checking that the number threads is square number
    int t = sqrt(threadCount);
    if (threadCount != t * t) {
        if (myRank == 0)
            printf("[Error] Number of OpenMP Threads must be square of an "
                   "integer.\n");
        MPI_Finalize();
        return 1;
    }
#endif

    // Getting side elements (N) of the input A and B matrices
    int N = atoi(argv[1]);

    // Calculating side elements per MPI rank
    int n = N / p;

    // Allocating local A,B,C submatrices and copies of A and B
    double* A = (double*)malloc(n * n * sizeof(double));
    double* B = (double*)malloc(n * n * sizeof(double));
    double* tmpA = (double*)malloc(n * n * sizeof(double));
    double* tmpB = (double*)malloc(n * n * sizeof(double));
    double* C = (double*)malloc(n * n * sizeof(double));

    /*****************************************************************
     * TODO: Calculate which position in the 2D processor corresponds
     * to the current MPI rank
     ****************************************************************/

    // Based on the grid, get my current X and Y positions
    int myRankY = 0;
    int myRankX = 0;

    /*****************************************************************
     * TODO: Calculate the MPI ranks with whom we will exchange
     * A and B submatrices.
     ****************************************************************/

    int rankSendB = 0;
    int rankRecvB = 0;
    int rankSendA = 0;
    int rankRecvA = 0;

    // Initializing values of A and B submatrices, initially shifted as per
    // Cannon's algorithm indication
    double *aptr, *bptr;
    aptr = A;
    bptr = B;

    for (int i = n * myRankX; i < n * (myRankX + 1); i++)
        for (int j = n * myRankY; j < n * (myRankY + 1); j++) {
            *aptr = 1.0 / ((i + myRankY * n) % N + j + 1);
            *bptr = 1.0 / ((j + myRankX * n) % N + i + 1);
            aptr++;
            bptr++;
        }

    // Initializing result submatrix (C) with all zeros
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i * n + j] = 0.0;

    // Starting to measure running time from this point
    MPI_Barrier(MPI_COMM_WORLD);
    double t_init = MPI_Wtime();

#ifdef HYBRID

    // If using MPI+OpenMP, multiply the A and B submatrices with cannonOMP
    ompCannon(A, B, C, n, t);

#else

    // If using only MPI, multiply the A and B submatrices with dgemm
    dgemm_("N", "N", &n, &n, &n, &one, A, &n, B, &n, &one, C, &n);

#endif

    for (int step = 1; step < p; step++) {
        /*****************************************************************
         * TODO: use MPI to exchange submatrices here.
         ****************************************************************/

        // Exchanging pointers between the received tmpA and tmpB submatrices
        // and the ones we operate with (A and B) when calling dgemm
        std::swap(A, tmpA);
        std::swap(B, tmpB);

#ifdef HYBRID

        // If using MPI+OpenMP, multiply the A and B submatrices with cannonOMP
        ompCannon(A, B, C, n, t);

#else

        // If using only MPI, multiply the A and B submatrices with dgemm
        dgemm_("N", "N", &n, &n, &n, &one, A, &n, B, &n, &one, C, &n);

#endif
    }

    // Barrier that waits for everyone to finish before taking the final timing
    MPI_Barrier(MPI_COMM_WORLD);
    double t_end = MPI_Wtime();

    // Verification phase
    int error = false;

#ifdef VERIFY
    if (myRank == 0)
        printf("[Info] Computation completed. Now verifying...\n");
    double tolerance = 1e-6;
    double* cptr = C;
    for (int i = n * myRankX; i < n * (myRankX + 1) && !error; i++)
        for (int j = n * myRankY; j < n * (myRankY + 1) && !error;
             j++, cptr++) {
            double tmp = 0;
            for (int k = 0; k < N; k++)
                tmp += 1.0 / ((i + k + 1) * (k + j + 1));
            error = fabs(*cptr - tmp) > tolerance;
        }
    int tempErr = error;
    MPI_Reduce(&tempErr, &error, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
#endif

    // Calculating performance
    double execTime = t_end - t_init;
    double gflops = ((2e-9) * N * N * N) / execTime;

    // Printing result
    if (myRank == 0) {
        if (error) {
            printf("[Error] Verification Failed!\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        printf("[Success] Execution time: %.3fs (GFlop/s: %.4f)\n", execTime,
               gflops);
    }

    free(A);
    free(B);
    free(tmpA);
    free(tmpB);
    free(C);

    MPI_Finalize();
    return 0;
}
