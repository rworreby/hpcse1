#include <vector>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>

std::vector<double> Ax_row(std::vector<double> &A, std::vector<double> &x){

    size_t N = x.size();
    std::vector<double> y(N);

    //TOTO: Question 2a: Matrix vector multiplication with a row major matrix


    return y;
}


std::vector<double> Ax_col(std::vector<double> &A, std::vector<double> &x){

    size_t N = x.size();
    std::vector<double> y(N);

    //TOTO: Question 2a: Matrix vector multiplication with a column major matrix

    return y;
}


double benchmark_Ax(std::vector<double> &A, std::vector<double> &x, bool row_major, double Ns){

    double times = 0;

    for( size_t i=0; i<Ns; i++){
        auto t1 = std::chrono::system_clock::now();
        //TOTO: Question 2a: Call the function to be benchmarked
        if( row_major==true ){

        }
        else{

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

    // TODO: Question 2a: Initialize matrices and vector
    //       store A as row major and B as column major
    

    std::cout << "Working with matrix of dimension " << N << std::endl;

    std::cout << "A*x (row major)." << std::endl;
    double times1 = benchmark_Ax(A,x,true,Ns);

    std::cout << "A*x (column major)." << std::endl;
    double times2 = benchmark_Ax(B,x,false,Ns);

    std::cout << "-----------------\n";
    std::cout << "Speedup " << times1/times2 << std::endl;


    return 0;
}
