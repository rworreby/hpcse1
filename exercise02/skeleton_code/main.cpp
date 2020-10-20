#include <omp.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>

// Integrand
inline double F(double x, double y)
{
    if (x * x + y * y < 1.) { // inside unit circle
        return 4.;
    }
    return 0.;
}

// Method 0: serial
double C0(size_t n)
{
    // random generator with seed 0
    std::default_random_engine g(0);
    // uniform distribution in [0, 1]
    std::uniform_real_distribution<double> u;

    double s = 0.; // sum
    for (size_t i = 0; i < n; ++i) {
        double x = u(g);
        double y = u(g);
        s += F(x, y);
    }
    return s / n;
}

// Method 1: openmp, no arrays
// TODO: Question 1a.1
double C1(size_t n)
{
    return 1.;
}

// Method 2, only `omp parallel for reduction`, arrays without padding
// TODO: Question 1a.2
double C2(size_t n)
{
    return 1.;
}

// Method 3, only `omp parallel for reduction`, arrays with padding
// TODO: Question 1a.3
double C3(size_t n)
{
    return 1.;
}

// Returns integral of F(x,y) over unit square (0 < x < 1, 0 < y < 1).
// n: number of samples
// m: method
double C(size_t n, size_t m)
{
    switch (m) {
    case 0:
        return C0(n);
    case 1:
        return C1(n);
    case 2:
        return C2(n);
    case 3:
        return C3(n);
    default:
        printf("Unknown method '%ld'\n", m);
        abort();
    }
}

int main(int argc, char* argv[])
{
    // default number of samples
    const size_t ndef = 1e8;

    if (argc < 2 || argc > 3 || std::string(argv[1]) == "-h") {
        fprintf(stderr, "usage: %s METHOD [N=%ld]\n", argv[0], ndef);
        fprintf(stderr, "Monte-Carlo integration with N samples.\n\
METHOD:\n\
0: serial\n\
1: openmp, no arrays\n\
2: `omp parallel for reduction`, arrays without padding\n\
3: `omp parallel for reduction`, arrays with padding\n");
        return 1;
    }

    // method
    size_t m = atoi(argv[1]);
    // number of samples
    size_t n = (argc > 2 ? atoi(argv[2]) : ndef);
    // reference solution
    double ref = 3.14159265358979323846;

    double wt0 = omp_get_wtime();
    double res = C(n, m);
    double wt1 = omp_get_wtime();

    printf("res:  %.20f\nref:  %.20f\nerror: %.20e\ntime: %.20f\n", res, ref,
           res - ref, wt1 - wt0);

    return 0;
}
