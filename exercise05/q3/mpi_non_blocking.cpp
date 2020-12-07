#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <mpi.h>

static double gammaInit(double x)
{
    return 4 * x / std::sqrt(1 - 4 * x * x);
}

static void initialConditions(double start, double end, std::vector<double>& x,
                              std::vector<double>& y,
                              std::vector<double>& gamma)
{
    const int n = (int)gamma.size();
    const double dx = (end - start) / n;

    for (int i = 0; i < n; ++i) {
        x[i] = start + dx * (i + 0.5);
        y[i] = 0;
        gamma[i] = dx * gammaInit(x[i]);
    }
}

static void computeVelocities(MPI_Comm comm, double epsSq,
                              const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& gamma,
                              std::vector<double>& u, std::vector<double>& v)
{
    // TODO: perform multi pass to compute interactions and update the local
    // velocities with non blocking communication
}

static void forwardEuler(double dt, const std::vector<double>& u,
                         const std::vector<double>& v, std::vector<double>& x,
                         std::vector<double>& y)
{
    for (int i = 0; i < (int)x.size(); ++i) {
        x[i] += dt * u[i];
        y[i] += dt * v[i];
    }
}

static void dumpToCsv(int step, const std::vector<double>& x,
                      const std::vector<double>& y,
                      const std::vector<double>& gamma)
{
    char fname[128];
    sprintf(fname, "config_%05d.csv", step);

    FILE* f = fopen(fname, "wb");
    fprintf(f, "x,y,gamma\n");
    for (int i = 0; i < (int)x.size(); ++i) {
        fprintf(f, "%g,%g,%g\n", x[i], y[i], gamma[i]);
    }
    fclose(f);
}

static void dumpToCsv(MPI_Comm comm, int step, const std::vector<double>& x,
                      const std::vector<double>& y,
                      const std::vector<double>& gamma)
{
    int rank, nranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    std::vector<double> xAll, yAll, gammaAll;

    // TODO Gather the data on rank zero before dumping to the csv files.

    if (rank == 0)
        dumpToCsv(step, xAll, yAll, gammaAll);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    if (argc != 2) {
        fprintf(stderr, "usage: %s <total number of particles>\n", argv[0]);
        exit(1);
    }

    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    const int nglobal = std::atoi(argv[1]);

    if (nglobal % nranks != 0) {
        fprintf(stderr,
                "expected n to be a multiple of the number of ranks.\n");
        exit(1);
    }

    // TODO initialize the data for each rank.

    const int n = nglobal;      // TODO
    const double extents = 1.0; // TODO

    const double startX = -0.5; // TODO
    const double endX = startX + extents;

    const double dt = 1e-4;
    const double epsSq = 1e-3;
    const double endTime = 2.5;
    const double dumpFreq = 0.1;

    const int dumpEvery = dumpFreq / dt;
    const int numSteps = endTime / dt;

    // state of the simulation
    std::vector<double> x(n);
    std::vector<double> y(n);
    std::vector<double> gamma(n);

    // workspace
    std::vector<double> u(n);
    std::vector<double> v(n);

    initialConditions(startX, endX, x, y, gamma);

    for (int step = 0; step < numSteps; ++step) {
        if (step % dumpEvery == 0) {
            const int id = step / dumpEvery;
            dumpToCsv(comm, id, x, y, gamma);
        }

        computeVelocities(comm, epsSq, x, y, gamma, u, v);
        forwardEuler(dt, u, v, x, y);
    }

    MPI_Finalize();

    return 0;
}
