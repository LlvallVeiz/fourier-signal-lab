/* fft.c -- Iterative radix-2 Cooley-Tukey FFT and naive DFT reference.
 *
 * The FFT follows the classic decimation-in-time formulation of Cooley
 * and Tukey (1965): an N-point DFT is split into the DFTs of the
 * even- and odd-indexed samples, recursed until length 1, then
 * recombined with twiddle factors. The iterative version below does
 * the same work in place: a bit-reversal permutation followed by
 * log2(N) butterfly stages, O(N log N) operations and O(1) extra
 * storage. The DFT reference is the direct O(N^2) definition.
 */

#include "fft.h"

#include <math.h>

int fft_is_pow2(size_t n) {
    return n != 0 && (n & (n - 1)) == 0;
}

size_t fft_next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}

/* Reorder x in place so index i moves to its bit-reversed position. */
static void bit_reverse(cx_t *x, size_t n) {
    size_t j = 0;
    for (size_t i = 1; i < n; i++) {
        size_t bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            cx_t tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
    }
}

int fft_forward(cx_t *x, size_t n) {
    if (x == NULL || !fft_is_pow2(n)) {
        return -1;
    }
    bit_reverse(x, n);
    for (size_t len = 2; len <= n; len <<= 1) {
        /* Twiddle factor for this stage: exp(-2*pi*i / len). */
        const double ang = -2.0 * FL_PI / (double)len;
        const cx_t wlen = cx_make(cos(ang), sin(ang));
        const size_t half = len / 2;
        for (size_t i = 0; i < n; i += len) {
            cx_t w = cx_make(1.0, 0.0);
            for (size_t j = 0; j < half; j++) {
                const cx_t u = x[i + j];
                const cx_t v = cx_mul(x[i + j + half], w);
                x[i + j] = cx_add(u, v);
                x[i + j + half] = cx_sub(u, v);
                w = cx_mul(w, wlen);
            }
        }
    }
    return 0;
}

int fft_inverse(cx_t *x, size_t n) {
    if (x == NULL || !fft_is_pow2(n)) {
        return -1;
    }
    /* Inverse via conjugation: IFFT(x) = conj(FFT(conj(x))) / N. */
    for (size_t i = 0; i < n; i++) {
        x[i] = cx_conj(x[i]);
    }
    if (fft_forward(x, n) != 0) {
        return -1;
    }
    const double scale = 1.0 / (double)n;
    for (size_t i = 0; i < n; i++) {
        x[i] = cx_scale(cx_conj(x[i]), scale);
    }
    return 0;
}

int dft_forward(const cx_t *in, cx_t *out, size_t n) {
    if (in == NULL || out == NULL || n == 0) {
        return -1;
    }
    for (size_t k = 0; k < n; k++) {
        double re = 0.0;
        double im = 0.0;
        for (size_t t = 0; t < n; t++) {
            const double ang = -2.0 * FL_PI * (double)k * (double)t / (double)n;
            const double c = cos(ang);
            const double s = sin(ang);
            re += in[t].re * c - in[t].im * s;
            im += in[t].re * s + in[t].im * c;
        }
        out[k] = cx_make(re, im);
    }
    return 0;
}

int dft_inverse(const cx_t *in, cx_t *out, size_t n) {
    if (in == NULL || out == NULL || n == 0) {
        return -1;
    }
    for (size_t t = 0; t < n; t++) {
        double re = 0.0;
        double im = 0.0;
        for (size_t k = 0; k < n; k++) {
            const double ang = 2.0 * FL_PI * (double)k * (double)t / (double)n;
            const double c = cos(ang);
            const double s = sin(ang);
            re += in[k].re * c - in[k].im * s;
            im += in[k].re * s + in[k].im * c;
        }
        out[t] = cx_make(re / (double)n, im / (double)n);
    }
    return 0;
}
