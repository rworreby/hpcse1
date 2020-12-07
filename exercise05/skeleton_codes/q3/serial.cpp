#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

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

static void computeVelocities(double epsSq, const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& gamma,
                              std::vector<double>& u, std::vector<double>& v)
{
    // TODO compute the interactions between the particles and write the result
    // in the velocity vectors u and v
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

int main(int argc, char** argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s <total number of particles>\n", argv[0]);
        exit(1);
    }

    const int n = std::atoi(argv[1]);

    const double startX = -0.5;
    const double endX = 0.5;

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
            dumpToCsv(id, x, y, gamma);
        }

        computeVelocities(epsSq, x, y, gamma, u, v);
        forwardEuler(dt, u, v, x, y);
    }

    return 0;
}
