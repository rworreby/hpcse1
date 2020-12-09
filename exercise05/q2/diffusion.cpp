#include <fstream>
#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <math.h>

struct Diagnostics {
    double time;
    double concentration;

    Diagnostics(double time, double concentration) : time(time), concentration(concentration) {}
};

struct Diffusion
{
    double D, L;  //diffusion constant and domain length
    int N;        //grid points per direction (whole grid is NxN)
    int local_N;  //number of rows of this process

    double h, dt;  //grid spacing and timestep
    double aux;   //auxiliary variable

    std::vector<double> c;    //solution vector
    std::vector<double> c_tmp;

    int rank, size; //MPI rank and total number of ranks

    std::vector<Diagnostics> diag;


    Diffusion(double D, double L, int N, int rank, int size) : D(D), L(L), N(N), rank(rank), size(size)
    {
        h = L / (N - 1);
        dt = h*h/(4.0*D); //this is the largest possible timestep (larger values lead to instabilities)

        local_N = N / size;
        if (rank == size - 1) local_N += N % size; //Correction for the last process

        c.resize((local_N+2)*(N+2), 0.0); //+2 for the ghost cells
        c_tmp.resize((local_N+2)*(N+2), 0.0);

        aux = dt * D / (h*h);
        initialize_density();
    }

    void advance()
    {
        int prev_rank{ rank - 1 };
        int next_rank{ rank + 1 };

        if (prev_rank < 0){
            prev_rank = MPI_PROC_NULL;
        }
        if (next_rank >= size){
            next_rank = MPI_PROC_NULL;
        }

        /* Signature of MPI_Sendrecv
        int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                        int dest, int sendtag,
                        void *recvbuf, int recvcount, MPI_Datatype recvtype,
                        int source, int recvtag,
                        MPI_Comm comm, MPI_Status *status)
        */
        MPI_Sendrecv(&c[(N+2)*local_N+1], N, MPI_DOUBLE,
                     next_rank, 0,
                     &c[1], N, MPI_DOUBLE,
                     prev_rank, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&c[(N+2)+1], N, MPI_DOUBLE,
                     prev_rank, 0,
                     &c[(N+2)*(local_N+1)+1], N, MPI_DOUBLE,
                     next_rank, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        /* Central differences in space, forward Euler in time, Dirichlet BCs */
        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j)
                c_tmp[i * (N + 2) + j] =
                    c[i * (N + 2) + j] +
                    aux * (c[i * (N + 2) + (j + 1)] + c[i * (N + 2) + (j - 1)] +
                           c[(i + 1) * (N + 2) + j] + c[(i - 1) * (N + 2) + j] -
                           4 * c[i * (N + 2) + j]);

        /* Use swap instead of c = c_tmp. This is much more efficient,
           because it does not copy element by element, just replaces storage
           pointers. */
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

        MPI_Reduce(rank ? &amount : MPI_IN_PLACE, &amount, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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

        /* Find max and min density values */
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

        MPI_Allreduce(MPI_IN_PLACE, &max_c, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &min_c, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        double epsilon = 1e-8;
        double dh = (max_c - min_c + epsilon) / M;

        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j) {
                int bin = (c[i * (N + 2) + j] - min_c) / dh;
                hist[bin]++;
            }

        std::vector<int> g_hist(M, 0);

        MPI_Reduce(hist.data(), g_hist.data(), g_hist.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

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

    void initialize_density()
    {
        int gi; // global index
        double bound = 0.25 * L;

        for (int i = 0; i < local_N; ++i) {
            gi = rank * (N / size) + i; // convert local index to global index

            for (int j = 0; j < N; ++j) {
                if (fabs(gi * h - 0.5*L) < bound && fabs(j * h - 0.5*L) < bound)
                    c[(i+1) * (N+2) + (j+1)] = 1;
                else
                    c[(i+1) * (N+2) + (j+1)] = 0;
            }
        }
    }

};

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " D L N \n";
        std::cerr << "D: diffusion constant\n";
        std::cerr << "L: and domain length\n";
        std::cerr << "N: grid points per direction (whole grid is NxN)" << '\n';
        return 1;
    }


    int rank,size;
    rank = 0;
    size = 1;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const double D = std::stod(argv[1]);
    const double L = std::stod(argv[2]);
    const int N = std::stoul(argv[3]);

    if (rank == 0)
        printf("Running Diffusion 2D on a %d x %d grid with %d ranks.\n", N, N, size);

    Diffusion system(D, L, N, rank, size);
    system.compute_diagnostics(0);
    for (int step = 0; step < 10000; ++step) {
        system.advance();
        system.compute_diagnostics(system.dt * step);
    }
    system.compute_histogram();
    if (rank == 0)
        system.write_diagnostics("diagnostics.dat");

    MPI_Finalize();

    return 0;
}
