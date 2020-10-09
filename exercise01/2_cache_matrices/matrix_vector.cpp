#include <vector>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>

std::vector<double> Ax_row(std::vector<double> &A, std::vector<double> &x){

    size_t N = x.size();
    std::vector<double> y(N, 0);

    for (size_t i = 0; i < N; i++) {
        size_t const i_N{ i*N };
        double tmp{ 0.0 };
        double x_i{ x[i] };
        for (size_t j = 0; j < N; j++) {
            tmp += A[i_N + j] * x_i;
        }
        y[i] = tmp;
    }

    return y;
}


std::vector<double> Ax_col(std::vector<double> &A, std::vector<double> &x){

    size_t N = x.size();
    std::vector<double> y(N, 0);

    for (size_t i = 0; i < N; i++) {
        size_t const i_N{ i*N };
        for (size_t j = 0; j < N; j++) {
            y[j] += A[j*N + i] * x[j];
        }
    }
    return y;
}


double benchmark_Ax(std::vector<double> &A, std::vector<double> &x, bool row_major, double Ns){
    double times = 0;

    for(size_t i=0; i<Ns; i++){
        auto t1 = std::chrono::system_clock::now();
        if(row_major==true){
            Ax_row(A, x);
        }
        else{
            Ax_col(A, x);
        }
        auto t2 = std::chrono::system_clock::now();
        times += std::chrono::duration<double>(t2-t1).count();
    }
    std::cout << "Done in total " << times << " -- average " << times/Ns << '\n';

    return times/Ns;
}


int main(int argc, char const *argv[])
{
    if(argc<3){
        std::cout << "Usage: " << argv[0]
                  << " [N|matrix dimension] [Ns|number of iterations]"
                  << std::endl;
        exit(0);
    }

    size_t N = atoi(argv[1]);
    size_t Ns = atoi(argv[2]);
    std::vector<double> A(N*N), B(N*N), x(N);

    for (size_t i = 0; i < N; ++i) {
        x[i] = i;
        for (size_t j = 0; j < N; ++j) {
            A[j*N + i] = i + j;
            B[i*N + j] = i + j;
        }
    }

    std::cout << "Working with matrix of dimension " << N << std::endl;

    std::cout << "A*x (row major)." << std::endl;
    double times1 = benchmark_Ax(A,x,true,Ns);

    std::cout << "A*x (column major)." << std::endl;
    double times2 = benchmark_Ax(B,x,false,Ns);

    std::cout << "-----------------\n";
    std::cout << "Speedup " << times1/times2 << std::endl;


    return 0;
}
