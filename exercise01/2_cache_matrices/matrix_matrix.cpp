#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <string>


void AB(std::vector<double> A, std::vector<double> B){

    size_t N = sqrt(A.size());
    std::vector<double> C(N*N);

    // TODO: Question 2c: Straightforward matrix-matrix multiplication

}


void AB_block_row(std::vector<double> A, std::vector<double> B, size_t blockSize){

    size_t N = sqrt(A.size());
    std::vector<double> C(N*N);

    // TODO: Question 2c: Block matrix-matrix multiplication - B in row major

}


void AB_block_col(std::vector<double> A, std::vector<double> B, size_t blockSize){

    size_t N = sqrt(A.size());
    std::vector<double> C(N*N);

    // TODO: Question 2c: Block matrix-matrix multiplication - B in column major

}


double benchmark_AB(std::vector<double> A, std::vector<double> B, size_t mode, size_t blockSize, size_t Ns){

    size_t N = sqrt(A.size());
    double times = 0;

    // TODO: Check that the matrix size is divided by the blockSize when mode==2 or 3
    if((mode==2 or mode==3) &&  N%blockSize!=0){
        std::cout << "Error: the size of the matrix " << N
                  << " should be divided by the blockSize variable "
                  << blockSize << std::endl;
        exit(1);
    }

    for(size_t i=0; i<Ns; i++){
        auto t1 = std::chrono::system_clock::now();
        // TODO: Question 2c: Call the function to be benchmarked
        if( mode==1 ){

        }
        else if( mode==2 ){

        }
        else if( mode==3 ){

        }
        auto t2 = std::chrono::system_clock::now();
        times += std::chrono::duration<double>(t2-t1).count();
    }
    std::cout << "Done in total " << times << " -- average " << times/Ns << '\n';

    return times/Ns;

}


int main(){
    std::vector<int> matrixSize{ 256, 512, 1024, 2048 };
    size_t M = matrixSize.size();

    std::vector<size_t> blockSize{ 2, 4, 8, 16, 32, 64, 128 };
    size_t Bs = blockSize.size();

    size_t Ns = 5;

    std::vector<double> times1(M);
    std::vector< std::vector<double>> times2(Bs, std::vector<double>(M) );
    std::vector< std::vector<double>> times3(Bs, std::vector<double>(M) );

    std::vector<double> A, B, C;


    for(size_t m=0; m<M; m++){
        std::cout << "Working with matrices of size " << matrixSize[m] << '\n';
        std::cout << "---------------------------------------------\n";

        size_t N = matrixSize[m];

        // TODO: Question 2c: Initialize matrices
        //       store A and B as row major and C as column major

        std::cout << "Start C=A*B (non optimized)." << '\n';
        times1[m] = benchmark_AB( A, B, 1, 0, Ns );

        std::cout << "---------------------------------------------\n";

        for(size_t b=0; b<Bs; b++){
            std::cout << "Start C=A*B (optimized, row major, block size="
                      << blockSize[b] << std::endl;
            times2[b][m] = benchmark_AB( A, B, 2, blockSize[b], Ns );
        }

        std::cout << "---------------------------------------------\n";

        for(size_t b=0; b<Bs; b++){
            std::cout << "Start C=A*B (optimized, column major, block size="
                      << blockSize[b] << std::endl;
            times3[b][m] = benchmark_AB( A, C, 3, blockSize[b], Ns );
        }

        std::cout << "==================================================\n";
    }


    FILE *fp=nullptr;
    fp = fopen("matrix_matrix_times.txt","w");
    // write header to the file
    std::string header = " N   unopt ";
    for(size_t b=0; b<Bs; b++)
    header = header + "  br_" + std::to_string(blockSize[b]);
    for(size_t b=0; b<Bs; b++)
    header = header + "  bc" + std::to_string(blockSize[b]);
    header = header + "\n";
    fprintf(fp,"%s",header.c_str());

    for(size_t m=0; m<M; m++){
        fprintf(fp,"%d %lf",matrixSize[m],times1[m]);
        for(size_t b=0; b<Bs; b++)
        fprintf(fp," %lf ",times2[b][m]);
        for(size_t b=0; b<Bs; b++)
        fprintf(fp," %lf ",times3[b][m]);
        fprintf(fp,"\n");
    }
    fclose(fp);


    return 0;
}
