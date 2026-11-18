/* test_fft.c -- Tests for the DFT reference and the radix-2 FFT.
 *
 * Hand-computed values used here:
 * - DFT of the impulse [1,0,0,0] is [1,1,1,1] (every bin sums the
 *   single 1 with a unit-magnitude twiddle).
 * - DFT of [0,1,0,0]: X[k] = exp(-2*pi*i*k/4), so X[1] = (0,-1).
 * - A sine of amplitude A at exactly bin k of an N-point DFT has
 *   |X[k]| = |X[N-k]| = A*N/2 and ~0 elsewhere.
 * - Parseval: sum(x^2) = (1/N) * sum(|X|^2).
 */

#include "test_framework.h"

#include <math.h>
#include <stdlib.h>

#include "complex.h"
#include "fft.h"

static double max_cx_diff(const cx_t *a, const cx_t *b, size_t n) {
    double m = 0.0;
    for (size_t i = 0; i < n; i++) {
        const double d = fmax(fabs(a[i].re - b[i].re), fabs(a[i].im - b[i].im));
        if (d > m) {
            m = d;
        }
    }
    return m;
}

/* Deterministic pseudo-random doubles in [-1, 1). */
static double next_rand(unsigned long *state) {
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    return (double)((*state >> 11) & 0x1FFFFFFFFFFFFFULL) /
               4503599627370496.0 * 2.0 -
           1.0;
}

void test_fft(void) {
    /* --- DFT of the impulse [1,0,0,0] is all ones. --- */
    {
        cx_t in[4] = {cx_make(1, 0), cx_make(0, 0), cx_make(0, 0),
                      cx_make(0, 0)};
        cx_t out[4];
        CHECK(dft_forward(in, out, 4) == 0);
        for (int k = 0; k < 4; k++) {
            CHECK_NEAR(out[k].re, 1.0, 1e-12);
            CHECK_NEAR(out[k].im, 0.0, 1e-12);
        }
    }

    /* --- DFT of [0,1,0,0]: X[k] = exp(-2*pi*i*k/4). --- */
    {
        cx_t in[4] = {cx_make(0, 0), cx_make(1, 0), cx_make(0, 0),
                      cx_make(0, 0)};
        cx_t out[4];
        CHECK(dft_forward(in, out, 4) == 0);
        CHECK_NEAR(out[0].re, 1.0, 1e-12);
        CHECK_NEAR(out[0].im, 0.0, 1e-12);
        CHECK_NEAR(out[1].re, 0.0, 1e-12);
        CHECK_NEAR(out[1].im, -1.0, 1e-12);
        CHECK_NEAR(out[2].re, -1.0, 1e-12);
        CHECK_NEAR(out[2].im, 0.0, 1e-12);
        CHECK_NEAR(out[3].re, 0.0, 1e-12);
        CHECK_NEAR(out[3].im, 1.0, 1e-12);
    }

    /* --- DFT inverse round-trips a hand value. --- */
    {
        cx_t in[4] = {cx_make(1, 0), cx_make(0, -1), cx_make(-1, 0),
                      cx_make(0, 1)};
        cx_t out[4];
        CHECK(dft_inverse(in, out, 4) == 0);
        CHECK_NEAR(out[0].re, 0.0, 1e-12);
        CHECK_NEAR(out[1].re, 1.0, 1e-12);
        CHECK_NEAR(out[2].re, 0.0, 1e-12);
        CHECK_NEAR(out[3].re, 0.0, 1e-12);
    }

    /* --- FFT agrees with the naive DFT (real input, N=64 and N=256). --- */
    {
        const size_t n = 64;
        cx_t in[64], ref[64], fast[64];
        unsigned long state = 12345;
        for (size_t i = 0; i < n; i++) {
            in[i] = cx_make(next_rand(&state), 0.0);
            fast[i] = in[i];
        }
        CHECK(dft_forward(in, ref, n) == 0);
        CHECK(fft_forward(fast, n) == 0);
        CHECK(max_cx_diff(ref, fast, n) < 1e-10);
    }
    {
        const size_t n = 256;
        cx_t in[256], ref[256], fast[256];
        unsigned long state = 999;
        for (size_t i = 0; i < n; i++) {
            /* Complex input exercises the full butterfly. */
            in[i] = cx_make(next_rand(&state), next_rand(&state));
            fast[i] = in[i];
        }
        CHECK(dft_forward(in, ref, n) == 0);
        CHECK(fft_forward(fast, n) == 0);
        CHECK(max_cx_diff(ref, fast, n) < 1e-9);
    }

    /* --- Pure sine exactly on a bin: |X[k]| = A*N/2. --- */
    {
        const size_t n = 64;
        const double amp = 2.0;
        cx_t x[64];
        for (size_t i = 0; i < n; i++) {
            x[i] = cx_make(amp * sin(2.0 * FL_PI * 8.0 * (double)i / 64.0),
                           0.0);
        }
        CHECK(fft_forward(x, n) == 0);
        CHECK_NEAR(cx_abs(x[8]), 64.0, 1e-9);  /* A*N/2 = 64 */
        CHECK_NEAR(cx_abs(x[56]), 64.0, 1e-9); /* mirror bin */
        for (size_t k = 0; k < n; k++) {
            if (k != 8 && k != 56) {
                CHECK(cx_abs(x[k]) < 1e-9);
            }
        }
    }

    /* --- FFT of zeros is zeros; FFT of a constant is DC only. --- */
    {
        cx_t x[8];
        for (int i = 0; i < 8; i++) {
            x[i] = cx_make(0.0, 0.0);
        }
        CHECK(fft_forward(x, 8) == 0);
        for (int i = 0; i < 8; i++) {
            CHECK_NEAR(cx_abs(x[i]), 0.0, 1e-15);
        }
        for (int i = 0; i < 8; i++) {
            x[i] = cx_make(3.0, 0.0);
        }
        CHECK(fft_forward(x, 8) == 0);
        CHECK_NEAR(x[0].re, 24.0, 1e-12);
        CHECK_NEAR(x[0].im, 0.0, 1e-12);
        for (int i = 1; i < 8; i++) {
            CHECK(cx_abs(x[i]) < 1e-12);
        }
    }

    /* --- IFFT round trip: FFT -> IFFT recovers the input. --- */
    {
        const size_t n = 128;
        cx_t x[128], orig[128];
        unsigned long state = 424242;
        for (size_t i = 0; i < n; i++) {
            orig[i] = cx_make(next_rand(&state), next_rand(&state));
            x[i] = orig[i];
        }
        CHECK(fft_forward(x, n) == 0);
        CHECK(fft_inverse(x, n) == 0);
        CHECK(max_cx_diff(orig, x, n) < 1e-12);
    }

    /* --- Parseval: sum(x^2) = (1/N) sum(|X|^2). --- */
    {
        const size_t n = 256;
        cx_t x[256];
        double time_sum = 0.0;
        unsigned long state = 777;
        for (size_t i = 0; i < n; i++) {
            const double v = next_rand(&state);
            x[i] = cx_make(v, 0.0);
            time_sum += v * v;
        }
        CHECK(fft_forward(x, n) == 0);
        double freq_sum = 0.0;
        for (size_t k = 0; k < n; k++) {
            freq_sum += x[k].re * x[k].re + x[k].im * x[k].im;
        }
        freq_sum /= (double)n;
        CHECK(fabs(time_sum - freq_sum) / time_sum < 1e-12);
    }

    /* --- Invalid sizes are rejected, not crashed on. --- */
    {
        cx_t x[8];
        CHECK(fft_forward(x, 100) == -1); /* not a power of two */
        CHECK(fft_forward(x, 0) == -1);
        CHECK(fft_forward(NULL, 8) == -1);
        CHECK(fft_inverse(x, 12) == -1);
        CHECK(dft_forward(NULL, x, 8) == -1);
        CHECK(dft_forward(x, NULL, 8) == -1);
        CHECK(dft_forward(x, x, 0) == -1);
    }

    /* --- Helpers: power-of-two predicates. --- */
    CHECK(fft_is_pow2(1) && fft_is_pow2(2) && fft_is_pow2(1024));
    CHECK(!fft_is_pow2(0) && !fft_is_pow2(3) && !fft_is_pow2(1000));
    CHECK(fft_next_pow2(1) == 1);
    CHECK(fft_next_pow2(1000) == 1024);
    CHECK(fft_next_pow2(1024) == 1024);
}
