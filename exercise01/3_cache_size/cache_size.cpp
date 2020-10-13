#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

const int kMinN = 1000 / sizeof(int); // From 1 KB
const int kMaxN = 20000000 / sizeof(int); // to 20 MB.
const int kNumSamples = 100;
const int kM = 100000000; // Operations per sample.
int a[kMaxN]; // Permutation array.


void sattolo(int *p, int N) {
    /*
    * Generate a random single-cycle permutation using Satollo's algorithm.
    *
    * https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle#Sattolo's_algorithm
    * https://danluu.com/sattolo/
    */
    for (int i = 0; i < N; ++i){
        p[i] = i;
    }
    for (int i = 0; i < N - 1; ++i){
        std::swap(p[i], p[i + 1 + rand() % (N - i - 1)]);
    }
}

double measure(int N, int mode) {
    std::chrono::time_point<std::chrono::steady_clock> t_start, t_end;
    if (mode == 0) {
        int *p = a;
        sattolo(p, N);

    }
    else if (mode == 1) {
        for (size_t i = 0; i < N; i++) {
            a[i] = i + 1;
        }

    }
    else if (mode == 2) {
        for (int i = 0; i < N; ++i) {
            (i + 16 < N) ? a[i] = i + 16 : a[i] = (i + 16) % N;
        }
    }

    // TODO: Question 1b: Traverse the list (make M jumps, starting from k = 0) and measure the execution time.
    t_start = std::chrono::steady_clock::now();

    unsigned int tmp = 0;
    for(unsigned int i = 0; i < kM; ++i) {
        tmp = a[tmp];
    }

    t_end = std::chrono::steady_clock::now();
    return static_cast<std::chrono::duration<double>>(t_end - t_start).count();
}

void run_mode(int mode) {
    /*
    * Run the measurement for many different values of N and output in a
    * format compatible with the plotting script.
    */
    printf("%9s  %9s  %7s  %7s\n", "N", "size[kB]", "t[s]", "op_per_sec[10^9]");
    for (size_t i = 0; i < kNumSamples; ++i) {
        // Generate N in a logarithmic scale.
        size_t N = static_cast<int>(kMinN *
            std::pow(
                static_cast<double>(kMaxN) / static_cast<double>(kMinN),
                static_cast<double>(i) / static_cast<double>(kNumSamples - 1)
            )
        );
        double t = measure(N, mode);
        printf("%9d  %9.1lf  %7.5lf  %7.6lf\n",
               N, N * sizeof(int) / 1024., t, kM / t * 1e-9
              );
        fflush(stdout);
    }
    printf("\n\n");
}

int main() {
    run_mode(1);
    // run_mode(1);
    // run_mode(2);

    return 0;
}
