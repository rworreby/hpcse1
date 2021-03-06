#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <string>


void AB(std::vector<double> A, std::vector<double> B){

    size_t N = sqrt(A.size());
    std::vector<double> C(N*N, 0);

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            for (size_t k = 0; k < N; k++) {
                C[i*N + j] += A[i*N + k] * B[k*N + j];
            }
        }
    }

    std::cout << "AB C[i]" << '\n';
    for (size_t i = 0; i < 8; i++) {
        std::cout << C[i] << ' ';
    }
    for (size_t i = 0; i < 8; i++) {
        std::cout << C[N*i] << ' ';
    }
    std::cout << '\n';
}


void AB_block_row(std::vector<double> A, std::vector<double> B, size_t block_size){

    size_t N = sqrt(A.size());
    std::vector<double> C(N*N, 0);

    for (size_t bi = 0; bi < N; bi+=block_size) {
        for (size_t bj = 0; bj < N; bj+=block_size) {
            for (size_t bk = 0; bk < N; bk+=block_size) {
                for (size_t i = 0; i < block_size; i++) {
                    for (size_t j = 0; j < block_size; j++) {
                        for (size_t k = 0; k < block_size; k++) {
                            C[(bi+i)*N + bj+j] += A[(bi+i)*N + bk+k] * B[(bk+k)*N + bj+j];
                        }
                    }
                }
            }
        }
    }

    std::cout << "AB_block_row C[i]" << '\n';
    for (size_t i = 0; i < 8; i++) {
        std::cout << C[i] << ' ';
    }
    for (size_t i = 0; i < 8; i++) {
        std::cout << C[N*i] << ' ';
    }
    std::cout << '\n';
}


void AB_block_col(std::vector<double> A, std::vector<double> B, size_t block_size){

    size_t N = sqrt(A.size());
    std::vector<double> C(N*N, 0);

    for (size_t bi = 0; bi < N; bi+=block_size) {
        for (size_t bj = 0; bj < N; bj+=block_size) {
            for (size_t bk = 0; bk < N; bk+=block_size) {
                for (size_t i = 0; i < block_size; i++) {
                    for (size_t j = 0; j < block_size; j++) {
                        for (size_t k = 0; k < block_size; k++) {
                            C[(bi+i)*N + bj+j] += A[(bi+i)*N + bk+k] * B[(bi+i)*N + bk+k];
                        }
                    }
                }
            }
        }
    }

    std::cout << "AB_block_col C[i]" << '\n';
    for (size_t i = 0; i < 8; i++) {
        std::cout << C[i] << ' ';
    }
    for (size_t i = 0; i < 8; i++) {
        std::cout << C[N*i] << ' ';
    }
    std::cout << '\n';
}


double benchmark_AB(std::vector<double> A, std::vector<double> B, size_t mode, size_t block_size, size_t Ns){

    size_t N = sqrt(A.size());
    double times = 0;

    if((mode==2 or mode==3) &&  N%block_size!=0){
        std::cout << "Error: the size of the matrix " << N
                  << " should be divided by the block_size variable "
                  << block_size << std::endl;
        exit(1);
    }

    for(size_t i=0; i<Ns; i++){
        auto t1 = std::chrono::system_clock::now();

        if( mode==1 ){
            AB(A, B);
        }
        else if( mode==2 ){
            AB_block_row(A, B, block_size);
        }
        else if( mode==3 ){
            AB_block_col(A, B, block_size);
        }

        auto t2 = std::chrono::system_clock::now();
        times += std::chrono::duration<double>(t2-t1).count();
    }
    std::cout << "Done in total " << times << " -- average " << times/Ns << '\n';

    return times/Ns;
}


int main(){
    std::vector<int> matrix_size{ 256, 512, 1024, 2048 };
    size_t M = matrix_size.size();

    std::vector<size_t> block_size{ 2, 4, 8, 16, 32, 64, 128 };
    size_t Bs = block_size.size();

    size_t Ns = 1; //5

    std::vector<double> times1(M);
    std::vector<std::vector<double>> times2(Bs, std::vector<double>(M));
    std::vector<std::vector<double>> times3(Bs, std::vector<double>(M));

    std::vector<double> A, B, C;


    for(size_t m=0; m<M; m++){
        std::cout << "Working with matrices of size " << matrix_size[m] << '\n';
        std::cout << "---------------------------------------------\n";

        size_t N = matrix_size[m];

        // TODO: Question 2c: Initialize matrices
        A.resize(N*N);
        B.resize(N*N);
        C.resize(N*N);

        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < N; ++j){
                A[i*N + j] = 2*i + j;
                B[i*N + j] = 2*i + j;
                C[j*N + i] = 2*i + j;
            }
        }

        std::cout << "Start C=A*B (non optimized)." << '\n';
        times1[m] = benchmark_AB( A, B, 1, 0, Ns );

        std::cout << "---------------------------------------------\n";

        for(size_t b=0; b<Bs; b++){
            std::cout << "Start C=A*B (optimized, row major, block size="
                      << block_size[b] << std::endl;
            times2[b][m] = benchmark_AB( A, B, 2, block_size[b], Ns );
        }

        std::cout << "---------------------------------------------\n";

        for(size_t b=0; b<Bs; b++){
            std::cout << "Start C=A*B (optimized, column major, block size="
                      << block_size[b] << std::endl;
            times3[b][m] = benchmark_AB( A, C, 3, block_size[b], Ns );
        }

        std::cout << "==================================================\n";
    }

    FILE *fp=nullptr;
    fp = fopen("matrix_matrix_times.txt","w");
    // write header to the file
    std::string header = "N,unopt";
    for(size_t b=0; b<Bs; b++)
    header = header + ",br_" + std::to_string(block_size[b]);
    for(size_t b=0; b<Bs; b++)
    header = header + ",bc" + std::to_string(block_size[b]);
    header = header + "\n";
    fprintf(fp,"%s",header.c_str());

    for(size_t m=0; m<M; m++){
        fprintf(fp,"%d,%lf",matrix_size[m],times1[m]);
        for(size_t b=0; b<Bs; b++){
            fprintf(fp,",%lf",times2[b][m]);
        }
        for(size_t b=0; b<Bs; b++){
            fprintf(fp,",%lf",times3[b][m]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    return 0;
}
