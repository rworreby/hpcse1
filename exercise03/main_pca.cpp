// #include "network/Network.h"
// #include "network/Optimizer.h"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <list>
#include <mkl_lapack.h>
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "code/utils.h"


int main (int argc, char** argv)
{
    // int D = 2;  // Data dimension
    // int N = 1024; // Number of training samples
    // int num_comp = 2; // Number of principal components
    // std::string data_name = "2D";
    // std::string scaler = "center";
    // std::string method_name = "PCA";
    // std::string data_path = "./data/" + data_name + "_dataset.txt";

    int D = 1850;  // Data dimension
    int N = 1280; // Number of training samples
    int num_comp = 10; // Number of principal components
    std::string data_name = "faces";
    std::string scaler = "standard";
    std::string method_name = "PCA";
    std::string data_path = "./data/"+data_name+"_dataset.txt";

    ///////////////////////////////////////////////////////////////////////////
    // Reading data. The data dimension is N x D.  The returned pointer points
    // to the data in row-major order. That is, if (i,j) corresponds to
    // the row and column index, respectively, you access the data with
    // data_{i,j} = data[i*D + j], where 0 <= i < N and 0 <= j < D.
    // The pointer points to a total memory allocation of N * D doubles.
    // In this example, i corresponds to the sample number n and j to the data
    // dimension (equivalently data[n*D+d]).
    ///////////////////////////////////////////////////////////////////////////

    std::cout << "Loading dataset from: " << data_path << std::endl;
    double* data = utils::loadDataset(data_path, N, D);

    /////////////////////////////////////////////////////////////////////////
    // PCA IMPLEMENTATION
    // 1. Transpose data (save them to new array)
    // 2. Compute mean and standard deviation of your data
    // 3. Normalize the data
    // 4. Build the covariance matrix
    // 5. Compute the eigenvalues and eigenvectors of the covariance matrix.
    //    Use LAPACK here.
    // 6. Extract the principal components and save them
    // 7. Reduce the dimensionality of the data by applying the compresion
    // 8. Report the compression ratio
    // 9. Reconstruct a reference image from its compressed form and save
    //    it in a .txt file
    /////////////////////////////////////////////////////////////////////////

    double time_start;
    double time_end;

    /////////////////////////////////////////////////////////////////////////
    //    1. Transpose data, i.e. convert it from row major to column major for
    //    easier future computation
    double* data_T = new double[N*D];

    utils::transposeData(data_T, data, N, D);

    /////////////////////////////////////////////////////////////////////////
    //    2. Compute mean and standard deviation of your data
    time_start = omp_get_wtime();
    double* data_mean = new double[D];
    double* data_std = new double[D];

    utils::computeMean(data_mean, data_T, N, D);
    utils::computeStd(data_std, data_mean, data_T, N, D);

    // Writing the mean and standard deviation to a file
    utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_mean.txt", data_mean, 1, D);
    utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_std.txt", data_std, 1, D);

    // Time tracking
    time_end = omp_get_wtime();
    std::cout << "MEAN/STD TIME = " << time_end-time_start << " seconds\n";

    /////////////////////////////////////////////////////////////////////////
    //    3. Normalize the data

    if(scaler.compare("standard") == 0){
        utils::standardizeColMajor(data_T, data_mean, data_std, N, D);
    }
    else if (scaler.compare("center") == 0){
        utils::centerDataColMajor(data_T, data_mean, N, D);
    }

    /////////////////////////////////////////////////////////////////////////
    //    4. Build the covariance matrix
    time_start = omp_get_wtime();
    double* data_cov = new double[D*D];

    utils::constructCovariance(data_cov, data_T, N, D);

    time_end = omp_get_wtime();
    std::cout << "COVARIANCE-MATRIX TIME = " << time_end-time_start << " seconds\n";

    /////////////////////////////////////////////////////////////////////////
    //    5. Compute the eigenvalues and eigenvectors of the covariance matrix.
    //    Use LAPACK here.
    //    Consult the interface given in the link
    //    http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html#ga442c43fca5493590f8f26cf42fed4044
    time_start = omp_get_wtime();

    char jobz = 'V'; // Compute both eigenvalues and orthonormal eigenvectors
    char uplo = 'U';
    int info, lwork;

    double* W = new double[D]; // Eigenvalues
    double* work = new double[2];

    // First call to dsyev_() with lwork = -1 to determine optimal size of lwork array (cheap call)
    lwork=-1;

    dsyev_(&jobz, &uplo, &D, data_cov, &D, W, work, &lwork, &info);

    // Allocate optimal workspace
    lwork = (int)work[0];
    delete[] work;
    work = new double[lwork];

    dsyev_(&jobz, &uplo, &D, data_cov, &D, W, work, &lwork, &info);

    // Upon completion data_cov contains the orthonormal eigenvectors of the covariance matrix stored in ROW major format: data_cov(j,k)=data_cov[j*D+k]
    // The eigenvalues are saved in ASCENDING order

    time_end = omp_get_wtime();
    std::cout << "DSYEV TIME = " << time_end-time_start << " seconds\n";
    // clean-up
    delete[] work;

    /////////////////////////////////////////////////////////////////////////
    //    6. Extract the principal components & eigenvalues and save them
    utils::reverseArray(W, D); // From ascending to descending order
    utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_eig.txt", W, 1, D);
    delete[] W;

    double* V = new double[num_comp * D];
    utils::getEigenvectors(V, data_cov, num_comp, D);
    delete[] data_cov;

    // Saving the PCA components
    utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_components.txt", V, num_comp, D);

    /////////////////////////////////////////////////////////////////////////
    //    7. Reduce the dimensionality of the data by applying the compresion
    time_start = omp_get_wtime();

    double* data_reduced = new double[N * num_comp]; // saved in normal row major order (not transpose!)
    utils::reduceDimensionality(data_reduced, V, data_T, N, D, num_comp);

    // SAVE THE DATA TO FILE
    utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_data_reduced.txt", data_reduced, num_comp, N);

    time_end = omp_get_wtime();
    std::cout << "PCREDUCED TIME = " << time_end-time_start << " seconds\n";


    /////////////////////////////////////////////////////////////////////////
    //    8. Report the compression ratio

    double dim_old{ static_cast<double>(D * N) };
    double dim_comp{ num_comp * D + num_comp * N + 2 * D };
    std::cout << "COMPRESSION RATIO = " << dim_old/dim_comp << std::endl;
    std::cout << "Floating point numbers used: " << '\n';
    std::cout << "Original dataset: " << dim_old << '\n';
    std::cout << "Compressed data: " << dim_comp << '\n';

    /////////////////////////////////////////////////////////////////////////
    //    9. Reconstruct the data from their compressed form and save them in a
    //    .txt file

    double* data_rec = new double[N * D];
    utils::reconstructDatasetRowMajor(data_rec, V, data_reduced, N, D, num_comp);

    /////////////////////////////////////////////////////////////////////////
    // Un-Normalize the data
    if(scaler.compare("standard") == 0){
        utils::inverseStandarizeDatasetRowMajor(data_rec, data_mean, data_std, N, D);
    }
    else if (scaler.compare("center") == 0){
        utils::inverseCenterDatasetRowMajor(data_rec, data_mean, N, D);
    }

    utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_data_reconstructed.txt", data_rec, N, D);

    // Clean-up
    delete[] data_rec;
    delete[] data_reduced;
    delete[] V;
    delete[] data_std;
    delete[] data_mean;
    delete[] data_T;
    delete[] data;
    return 0;
}
