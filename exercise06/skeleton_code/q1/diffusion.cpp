#include <fstream>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <vector>

struct Diagnostics {
    double time;
    double concentration;

    Diagnostics(double time, double concentration)
        : time(time), concentration(concentration)
    {
    }
};

struct Diffusion {
    double D, L; // diffusion constant and domain length
    int N;       // grid points per direction (whole grid is NxN)
    int local_N; // number of rows of this process

    double h, dt; // grid spacing and timestep
    double aux;   // auxiliary variable

    std::vector<double> c; // solution vector
    std::vector<double> c_tmp;

    int rank, size; // MPI rank and total number of ranks

    std::vector<Diagnostics> diag;

    Diffusion(double D, double L, int N, int rank, int size)
        : D(D), L(L), N(N), rank(rank), size(size)
    {
        h = L / (N - 1);
        dt = h * h / (4.0 * D); // this is the largest possible timestep (larger
                                // values lead to instabilities)

        local_N = N / size;
        if (rank == size - 1)
            local_N += N % size; // Correction for the last process

        c.resize((local_N + 2) * (N + 2), 0.0); //+2 for the ghost cells
        c_tmp.resize((local_N + 2) * (N + 2), 0.0);

        aux = dt * D / (h * h);
        initialize();
    }

    void advance()
    {
        // *** start MPI part ***
        MPI_Status status[2];

        int prev_rank = rank - 1;
        int next_rank = rank + 1;

        if (prev_rank < 0)
            prev_rank = MPI_PROC_NULL;
        if (next_rank >= size)
            next_rank = MPI_PROC_NULL;

        /* Exchange ALL necessary ghost cells with neighboring ranks */
        MPI_Sendrecv(&c[1 * (N + 2) + 1], N, MPI_DOUBLE, prev_rank, 0,
                     &c[(local_N + 1) * (N + 2) + 1], N, MPI_DOUBLE, next_rank,
                     0, MPI_COMM_WORLD, &status[0]);

        MPI_Sendrecv(&c[local_N * (N + 2) + 1], N, MPI_DOUBLE, next_rank, 1,
                     &c[0 * (N + 2) + 1], N, MPI_DOUBLE, prev_rank, 1,
                     MPI_COMM_WORLD, &status[1]);
        // *** end MPI part ***

        /* Central differences in space, forward Euler in time, Dirichlet BCs */
        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j)
                c_tmp[i * (N + 2) + j] =
                    c[i * (N + 2) + j] +
                    aux * (c[i * (N + 2) + (j + 1)] + c[i * (N + 2) + (j - 1)] +
                           c[(i + 1) * (N + 2) + j] + c[(i - 1) * (N + 2) + j] -
                           4 * c[i * (N + 2) + j]);
        using std::swap;
        swap(c_tmp, c);
    }

    void compute_diagnostics(const double t)
    {
        double amount = 0.0;

        /* Integration to compute total concentration */
        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j)
                amount += c[i * (N + 2) + j] * h * h;

        // *** start MPI part ***
        MPI_Reduce(rank == 0 ? MPI_IN_PLACE : &amount, &amount, 1, MPI_DOUBLE,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        // *** end MPI part ***

        if (rank == 0) {
            std::cout << "t = " << t << " amount = " << amount << '\n';
            diag.push_back(Diagnostics(t, amount));
        }
    }

    void write_diagnostics(const std::string& filename) const
    {
        std::ofstream out_file(filename, std::ios::out);
        for (const Diagnostics& d : diag)
            out_file << d.time << ' ' << d.concentration << '\n';
        out_file.close();
    }

    void compute_histogram()
    {
        /* Number of histogram bins */
        const int M = 10;
        std::vector<int> hist(M, 0);

        /* Find max and min concentration values */
        double max_c, min_c, c0;
        max_c = c[1 * (N + 2) + 1];
        min_c = c[1 * (N + 2) + 1];

        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j) {
                c0 = c[i * (N + 2) + j];
                if (c0 > max_c)
                    max_c = c0;
                if (c0 < min_c)
                    min_c = c0;
            }

        double lmin_c = min_c;
        double lmax_c = max_c;

        // *** start MPI part ***
        MPI_Allreduce(&lmin_c, &min_c, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&lmax_c, &max_c, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        // *** end MPI part ***

        double epsilon = 1e-8;
        double dh = (max_c - min_c + epsilon) / M;

        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j) {
                int bin = (c[i * (N + 2) + j] - min_c) / dh;
                hist[bin]++;
            }

        std::vector<int> g_hist(M, 0);
        // *** start MPI part ***
        MPI_Reduce(hist.data(), g_hist.data(), M, MPI_INT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        // *** end MPI part ***

        if (rank == 0) {
            printf("=====================================\n");
            printf("Output of compute_histogram():\n");
            int gl = 0;
            for (int i = 0; i < M; i++) {
                printf("bin[%d] = %d\n", i, g_hist[i]);
                gl += g_hist[i];
            }
            printf("Total elements = %d\n", gl);
        }

    } // end public

    void initialize()
    {
        int gi; // global index
        double bound = 0.25 * L;

        for (int i = 0; i < local_N; ++i) {
            gi = rank * (N / size) + i; // convert local index to global index

            for (int j = 0; j < N; ++j) {
                if (fabs(gi * h - 0.5 * L) < bound &&
                    fabs(j * h - 0.5 * L) < bound)
                    c[(i + 1) * (N + 2) + (j + 1)] = 1;
                else
                    c[(i + 1) * (N + 2) + (j + 1)] = 0;
            }
        }
    }

    void write_mpi_sequential(const std::string& filename) const
    {
        // Calculating elements per rank
        const int fullN = N + 2; // My part of the grid plus 2 ghost cells
        const int nlocal = (local_N) * (fullN);
        const size_t elementCount = size * nlocal;

        // Reserving incoming buffer for root rank
        std::vector<double> outputData;

        if (rank == 0)
            outputData.resize(elementCount);

        // Gathering all the data into the root rank
        MPI_Gather(&c[1 * fullN], nlocal, MPI_DOUBLE, outputData.data(), nlocal,
                   MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Saving data to file
        if (rank == 0) {
            auto file =
                fopen(filename.c_str(), "wb"); // Write-only, binary mode
            fwrite(outputData.data(), sizeof(outputData[0]), outputData.size(),
                   file);
            fclose(file);
        }
    }

    void write_mpi_parallel(const std::string& filename) const
    {
        // TODO: add your MPI I/O code here
    }
};

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " D L N \n";
        return 1;
    }

    int rank, size;
    // *** start MPI part ***
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // *** end MPI part ***

    const double D = std::stod(argv[1]);
    const double L = std::stod(argv[2]);
    const int N = std::stoul(argv[3]);

    if (rank == 0)
        printf("Running Diffusion 2D on a %d x %d grid with %d ranks.\n", N, N,
               size);

    Diffusion system(D, L, N, rank, size);

    double t_start = MPI_Wtime();
    for (int step = 0; step < 10; step++) {
        system.advance();
    }
    double t_end = MPI_Wtime();

    if (rank == 0)
        printf("[Success] Time taken: %.3fs.\n", t_end - t_start);

    system.compute_histogram();

    auto t0_start = MPI_Wtime();
    system.write_mpi_sequential("dump_sequential.dat");
    auto t0_end = MPI_Wtime();

    auto t1_start = MPI_Wtime();
    system.write_mpi_parallel("dump_parallel.dat");
    auto t1_end = MPI_Wtime();

    if (rank == 0) {
        printf("Time taken by sequential I/O: %.3fs.\n", t0_end - t0_start);
        printf("Time taken by parallel MPI I/O: %.3fs.\n", t1_end - t1_start);
    }

    // *** start MPI part ***
    MPI_Finalize();
    // *** end MPI part ***
    return 0;
}
