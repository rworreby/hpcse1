#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Grid in domain [0, 1] x [0, 1].
struct Grid {
    Grid(int N) : N(N), N_ext(N + 2), dx(1. / N), field_size(N_ext * N_ext) {}

    // Returns flat index from components.
    // ix,iy: indices from [-1, N], indices -1 and N refer to ghost cells
    size_t operator()(int ix, int iy) const
    {
        assert(-1 <= ix && ix <= N);
        assert(-1 <= iy && iy <= N);
        return (N + 2) * (iy + 1) + ix + 1;
    }

    const int N;             // cells per dimension
    const int N_ext;         // cells per dimension, including ghost cells
    const double dx;         // grid step
    const size_t field_size; // size of field data
};

using Field = std::vector<double>;

double getMax(const Field& field, const Grid& grid)
{
    double max = -std::numeric_limits<double>::max();
    for (int iy = 0; iy < grid.N; ++iy) {
        for (int ix = 0; ix < grid.N; ++ix) {
            max = std::max(max, field[grid(ix, iy)]);
        }
    }
    return max;
}

void initField(Field& field, const Grid& grid)
{
    field.assign(grid.field_size, 0);
    for (int iy = 0; iy < grid.N; ++iy) {
        for (int ix = 0; ix < grid.N; ++ix) {
            const double x = (ix + 0.5) * grid.dx;
            const double y = (iy + 0.5) * grid.dx;
            const double kx = 4 * M_PI;
            const double ky = 2 * M_PI;
            field[grid(ix, iy)] =
                std::pow(std::sin(x * kx) * std::sin(y * ky), 2);
        }
    }
}

void dumpField(std::string path, const Field& field, const Grid& grid, double t,
               int step)
{
    std::ofstream out(path);
    out << "# time: " << t << '\n';
    out << "# step: " << step << '\n';
    out << "# dim: " << grid.N << '\n';
    for (int iy = 0; iy < grid.N; ++iy) {
        for (int ix = 0; ix < grid.N; ++ix) {
            out << field[grid(ix, iy)] << ' ';
        }
        out << '\n';
    }
}

std::string dumpField(int idump, const Field& field, const Grid& grid, double t,
                      int step)
{
    char path[255];
    sprintf(path, "dump_%04d.dat", idump);
    dumpField(path, field, grid, t, step);
    return path;
}

class SolverExplicit
{
  public:
    SolverExplicit(const Grid& grid) : grid(grid), rhs_(grid.field_size) {}

    void advance(Field& phi, double dt, double D)
    {
        const double k = dt * D / (grid.dx * grid.dx);
        for (int iy = 0; iy < grid.N; ++iy) {
            for (int ix = 0; ix < grid.N; ++ix) {
                rhs_[grid(ix, iy)] =
                    k * (phi[grid(ix + 1, iy)] + phi[grid(ix - 1, iy)] +
                         phi[grid(ix, iy + 1)] + phi[grid(ix, iy - 1)] -
                         4 * phi[grid(ix, iy)]);
            }
        }
        for (int iy = 0; iy < grid.N; ++iy) {
            for (int ix = 0; ix < grid.N; ++ix) {
                phi[grid(ix, iy)] += rhs_[grid(ix, iy)];
            }
        }
    }

    const Grid grid;

  private:
    Field rhs_;
};

class SolverAdi
{
  public:
    SolverAdi(const Grid& grid)
        : grid(grid)
    {
    }

    // TODO: Implement the ADI scheme with the Thomas algorithm.
    // Check SolverExplicit to see how to loop over cells and access their
    // neighbors. You can define additional functions and member variables for
    // temporary buffers.

    void advance(Field& phi, double dt, double D)
    {
    }

    const Grid grid;

  private:
    Field rhs_;
    std::vector<double> tmp_;
};

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " N dt T [ndumps=2]\n";
        return 1;
    }

    const int N = std::stoul(argv[1]);                    // cells per dimension
    const double dt = std::stod(argv[2]);                 // timestep
    const double T = std::stod(argv[3]);                  // final time
    const int ndumps = argc > 4 ? std::stod(argv[4]) : 2; // number of dumps
    const double D = 1;                   // diffusion coefficient
    const int nsteps = std::ceil(T / dt); // number of time steps

    Grid grid(N);
    Field phi;
    initField(phi, grid);

    int idump = 0;

    std::ofstream fstat("stat.dat");
    fstat << "t max\n";

    SolverAdi solver(grid);
    for (int step = 0; step <= nsteps; ++step) {
        const double t = dt * step;

        if (step) {
            solver.advance(phi, dt, D);
        }

        fstat << t << ' ' << getMax(phi, grid) << '\n';
        if ((ndumps > 0 && step == 0) ||
            (ndumps > 1 && idump * nsteps / (ndumps - 1) == step)) {
            auto path = dumpField(idump++, phi, grid, t, step);
            printf("t=%f step=%d %s\n", t, step, path.c_str());
        }
    }
}
