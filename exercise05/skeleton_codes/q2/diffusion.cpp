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

    double h,dt;  //grid spacing and timestep
    double aux;   //auxiliary variable

    std::vector<double> c;    //solution vector
    std::vector<double> c_tmp;

    int rank,size; //MPI rank and total number of ranks

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

        // TODO: Implement Blocking MPI communication to exchange the ghost
        // cells on a periodic domain required to compute the central finite
        // differences below.

        // *** start MPI part ***

        // *** end MPI part ***

        /* Central differences in space, forward Euler in time, Dirichlet BCs */
        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j)
                c_tmp[i * (N + 2) + j] =
                    c[i * (N + 2) + j] +
                    aux * (c[i * (N + 2) + (j + 1)] + c[i * (N + 2) + (j - 1)] +
                           c[(i + 1) * (N + 2) + j] + c[(i - 1) * (N + 2) + j] -
                           4 * c[i * (N + 2) + j]);

        /* Use swap instead of rho_ = rho_tmp__. This is much more efficient,
           because it does not copy element by element, just replaces storage
           pointers. */
        using std::swap;
        swap(c_tmp, c);
    }

    void compute_diagnostics(const double t)
    {
        double ammount = 0.0;

        /* Integration to compute total concentration */
        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j)
                ammount += c[i * (N + 2) + j] * h * h;

        // TODO: Sum total ammount from all ranks

        // *** start MPI part ***

        // *** end MPI part ***

        if (rank == 0) {
            std::cout << "t = " << t << " ammount = " << ammount << '\n';
            diag.push_back(Diagnostics(t, ammount));
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

        // TODO: Compute the global min and max concentration values on this myRank and
        // store the result in min_rho and max_rho, respectively.

        // *** start MPI part ***

        // *** end MPI part ***

        double epsilon = 1e-8;
        double dh = (max_c - min_c + epsilon) / M;

        for (int i = 1; i <= local_N; ++i)
            for (int j = 1; j <= N; ++j) {
                int bin = (c[i * (N + 2) + j] - min_c) / dh;
                hist[bin]++;
            }

        // TODO: Compute the sum of the histogram bins over all ranks and store
        // the result in the array g_hist.  Only myRank 0 must print the result.
        std::vector<int> g_hist(M, 0);
        // *** start MPI part ***

        // *** end MPI part ***

        if (rank == 0) {
            printf("=====================================\n");
            printf("Output of compute_histogram():\n");
            int gl = 0;
            for (int i = 0; i < M; i++) {
                printf("bin[%d] = %d\n", i, hist[i]);
                gl += hist[i];
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
        return 1;
    }

    // TODO: Start-up the MPI environment and determine this process' myRank ID
    // as well as the total number of processes (=ranks) involved in the
    // communicator

    int rank,size;
    rank = 0;
    size = 1;
    // *** start MPI part ***




    // *** end MPI part ***

    const double D = std::stod(argv[1]);
    const double L = std::stod(argv[2]);
    const int N = std::stoul(argv[3]);

    if (rank == 0)
        printf("Running Diffusion 2D on a %d x %d grid with %d ranks.\n",N,N,size);

    Diffusion system(D, L, N, rank, size);
    system.compute_diagnostics(0);
    for (int step = 0; step < 10000; ++step) {
        system.advance();
        system.compute_diagnostics(system.dt * step);
    }
    system.compute_histogram();
    if (rank == 0)
        system.write_diagnostics("diagnostics.dat");
    // TODO: Shutdown the MPI environment
    // *** start MPI part ***

    // *** end MPI part ***
    return 0;
}
