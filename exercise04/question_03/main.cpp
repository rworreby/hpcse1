// Copyright 2020 ETH Zurich. All Rights Reserved.

#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <mpi.h>

inline long exact(const long N){
    return N * (N + 1) / 2;
}

void reduce_mpi(const int rank, long& sum){
    MPI_Reduce(rank ? &sum : MPI_IN_PLACE, &sum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
}

// PRE: size is a power of 2 for simplicity
void reduce_manual(int rank, int size, long& sum){
    // TODO f): Implement a tree based reduction using blocking point-to-point communication.
    for (int i = size/2; i > 0; i /= 2) {
        long int other{ 0 };
        if (rank < i){
            MPI_Recv(&other, 1, MPI_LONG, rank+i, rank+i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += other;
        }
        else if (rank <= i * 2){
            MPI_Send(&sum, 1, MPI_LONG, rank-i, rank, MPI_COMM_WORLD);
        }
    }
}


int main(int argc, char** argv){
    const long N = 1000000000;

    // Initialize MPI
    int rank, size;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // -------------------------
    // Perform the local sum:
    // -------------------------
    long sum = 0;

    // Determine work load per rank
    long N_per_rank = N / size;

    long N_start = rank * N_per_rank + 1;
    long N_end = N_start + N_per_rank - 1;

    // N_start + (N_start+1) + ... + (N_start+N_per_rank-1)
    for(long i = N_start; i <= N_end; ++i){
        sum += i;
    }

    // -------------------------
    // Reduction
    // -------------------------
    //reduce_mpi(rank, sum);
    reduce_manual(rank, size, sum);

    // -------------------------
    // Print the result
    // -------------------------
    if(rank == 0){
        std::cout << '\n';
        std::cout << std::left << std::setw(25) << "Final result (exact): " << exact(N) << std::endl;
        std::cout << std::left << std::setw(25) << "Final result (MPI): " << sum << std::endl;
    }

    MPI_Finalize();
    return 0;
}
