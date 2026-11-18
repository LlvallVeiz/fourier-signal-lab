/* fft.h -- Discrete Fourier transform: naive DFT reference and iterative
 * radix-2 Cooley-Tukey FFT, plus the inverse transform.
 *
 * Convention (matches numpy.fft):
 *   X[k] = sum_{n=0}^{N-1} x[n] * exp(-2*pi*i*k*n / N)   (forward)
 *   x[n] = (1/N) * sum_{k=0}^{N-1} X[k] * exp(+2*pi*i*k*n / N)
 *
 * The FFT is in-place; the DFT reference writes to a separate buffer.
 * All functions return 0 on success, -1 on invalid input (NULL pointer,
 * N == 0, or N not a power of two for the FFT).
 */

#ifndef FOURIER_LAB_FFT_H
#define FOURIER_LAB_FFT_H

#include <stddef.h>

#include "complex.h"

/* True (non-zero) when n is a power of two. */
int fft_is_pow2(size_t n);

/* In-place iterative radix-2 decimation-in-time FFT (Cooley-Tukey).
 * Requires n to be a power of two. */
int fft_forward(cx_t *x, size_t n);

/* In-place inverse FFT, including the 1/N scaling. */
int fft_inverse(cx_t *x, size_t n);

/* Naive O(N^2) DFT, any N >= 1. Reference implementation used to
 * validate the FFT and as a teaching baseline. */
int dft_forward(const cx_t *in, cx_t *out, size_t n);

/* Naive O(N^2) inverse DFT, including the 1/N scaling. */
int dft_inverse(const cx_t *in, cx_t *out, size_t n);

/* Smallest power of two >= n (n >= 1). */
size_t fft_next_pow2(size_t n);

#endif /* FOURIER_LAB_FFT_H */
