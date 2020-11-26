// Copyright 2020 ETH Zurich. All Rights Reserved.

#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <mpi.h>

inline long exact(const long N){
    // TODO b): Implement the analytical solution.
    return 0;
}

void reduce_mpi(const int rank, long& sum){
    // TODO e): Perform the reduction using blocking collectives.
}

// PRE: size is a power of 2 for simplicity
void reduce_manual(int rank, int size, long& sum){
    // TODO f): Implement a tree based reduction using blocking point-to-point communication.
}


int main(int argc, char** argv){
    const long N = 1000000;

    // Initialize MPI
    int rank, size;
    // TODO c): Initialize MPI and obtain the rank and the number of processes (size)

    // -------------------------
    // Perform the local sum:
    // -------------------------
    long sum = 0;

    // Determine work load per rank
    long N_per_rank = N / size;

    // TODO d): Determine the range of the subsum that should be calculated by this rank.
    long N_start;
    long N_end;

    // N_start + (N_start+1) + ... + (N_start+N_per_rank-1)
    for(long i = N_start; i <= N_end; ++i){
        sum += i;
    }

    // -------------------------
    // Reduction
    // -------------------------
    reduce_mpi(rank, sum);
    //reduce_manual(rank, size, sum);

    // -------------------------
    // Print the result
    // -------------------------
    if(rank == 0){
        std::cout << std::left << std::setw(25) << "Final result (exact): " << exact(N) << std::endl;
        std::cout << std::left << std::setw(25) << "Final result (MPI): " << sum << std::endl;
    }
    // Finalize MPI
    // TODO c): Finalize MPI

    return 0;
}
