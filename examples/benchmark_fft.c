/* benchmark_fft.c -- Throwaway benchmark (not part of the library):
 * times the FFT against the naive DFT to make the O(N log N) vs
 * O(N^2) gap concrete. Compile with:
 *   gcc -std=c11 -O2 -Isrc examples/benchmark_fft.c src/fft.c \
 *       src/signal.c -o /tmp/benchmark_fft -lm
 * or simply `make bench`.
 */

#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "complex.h"
#include "fft.h"
#include "signal.h"

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(void) {
    const size_t big = 1u << 20; /* 1,048,576 samples */
    double *x = (double *)malloc(big * sizeof(double));
    cx_t *spec = (cx_t *)malloc(big * sizeof(cx_t));
    if (x == NULL || spec == NULL) {
        fprintf(stderr, "benchmark: out of memory\n");
        return 1;
    }
    signal_noise(x, big, 1.0, 7);
    for (size_t i = 0; i < big; i++) {
        spec[i] = cx_make(x[i], 0.0);
    }
    double t0 = now_s();
    fft_forward(spec, big);
    const double fft_big = now_s() - t0;

    const size_t m = 1024;
    cx_t *in = (cx_t *)malloc(m * sizeof(cx_t));
    cx_t *out = (cx_t *)malloc(m * sizeof(cx_t));
    cx_t *fast = (cx_t *)malloc(m * sizeof(cx_t));
    if (in == NULL || out == NULL || fast == NULL) {
        fprintf(stderr, "benchmark: out of memory\n");
        return 1;
    }
    for (size_t i = 0; i < m; i++) {
        in[i] = cx_make(x[i], 0.0);
        fast[i] = in[i];
    }
    t0 = now_s();
    dft_forward(in, out, m);
    const double dft_small = now_s() - t0;
    t0 = now_s();
    fft_forward(fast, m);
    const double fft_small = now_s() - t0;

    printf("FFT  N=2^20 (1,048,576): %.3f s\n", fft_big);
    printf("DFT  N=1024:             %.6f s\n", dft_small);
    printf("FFT  N=1024:             %.6f s\n", fft_small);
    printf("FFT/DFT speed-up at N=1024: %.0fx\n", dft_small / fft_small);
    printf("(A DFT at N=2^20 would need ~(2^20/1024)^2 = 1,048,576x the\n"
           " N=1024 DFT time -- about 9 hours at this rate -- while the\n"
           " FFT does it in %.3f s.)\n",
           fft_big);

    free(x);
    free(spec);
    free(in);
    free(out);
    free(fast);
    return 0;
}
