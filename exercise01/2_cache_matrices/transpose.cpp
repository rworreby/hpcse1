#include <iostream>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <string>


void transpose(std::vector<double> A){

    size_t N = sqrt(A.size());
    std::vector<double> AT(N*N);
    // TODO: Question 2b: Straightforward matrix transposition

}


void transpose_block( std::vector<double> A, size_t block_size ){

    size_t N = sqrt(A.size());
    std::vector<double> AT(N*N);
    // TODO: Question 2b: Block matrix transposition

}


double benchmark_transpose(std::vector<double> A, size_t mode, size_t block_size, size_t Ns){

    size_t N = sqrt(A.size());
    double times = 0;

    // TODO: Check that the matrix size is divided by the block_size when mode==2
    if(mode==2 &&  N%block_size!=0){
        std::cout << "Error: the size of the matrix " << N
        << " should be divided by the block_size variable "
        << block_size << ".\n";
        exit(1);
    }

    for( size_t i=0; i<Ns; i++){
        auto t1 = std::chrono::system_clock::now();
        // TODO: Question 2b: Call the function to be benchmarked
        if(mode==1){

        }
        else if(mode==2){

        }
        auto t2 = std::chrono::system_clock::now();
        times += std::chrono::duration<double>(t2-t1).count();
    }
    std::cout << "Done in total " << times << " --  average "
              << times/Ns << std::endl;

    return times/Ns;

}


int main( )
{
    std::vector<int> matrix_size{ 1024, 2048, 4096 };
    size_t M = matrix_size.size();

    std::vector<size_t> block_size{ 2, 4, 8, 16, 32, 64, 128 };
    size_t B = block_size.size();

    size_t Ns = 2;

    std::vector<double> times1(M);
    std::vector< std::vector<double>> times2(B, std::vector<double>(M) );

    std::vector<double> A;

    // loop over matrix sizes
    for( size_t m=0; m<M; m++){

        std::cout << "Working with a matrix of size " << matrix_size[m] << ".\n";

        size_t N = matrix_size[m];

        // TODO: Initialize the matrix


        std::cout << "Start transposing (non optimized).\n";
        times1[m] = benchmark_transpose( A, 1, 0, Ns );

        // loop over block sizes
        for(size_t b=0; b<B; b++){

            std::cout << "Start transposing (optimized, block size="
                      << block_size[b] << ".\n";
            times2[b][m] = benchmark_transpose(A, 2, block_size[b], Ns);
        }

        std::cout << "==================================================\n";
    }


    // write results to a file
    FILE *fp=nullptr;
    fp = fopen("transpose_times.txt","w");
    // write header to the file
    std::string header = "# N   time_unoptimized ";
    for(size_t b=0; b<B; b++){
        header = header + "  block_" + std::to_string(block_size[b]);
    }
    header = header + "\n";
    fprintf(fp, "%s", header.c_str());
    for(size_t m=0; m<M; m++){
        fprintf(fp, "%d %lf", matrix_size[m], times1[m]);
        for(size_t b=0; b<B; b++){
            fprintf(fp, " %lf ", times2[b][m]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    return 0;
}
