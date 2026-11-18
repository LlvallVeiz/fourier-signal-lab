/* complex.h -- Minimal complex arithmetic for fourier-signal-lab.
 *
 * A tiny value type (no C99 <complex.h> dependency) so the FFT core
 * stays explicit and portable. All operations are static inline.
 */

#ifndef FOURIER_LAB_COMPLEX_H
#define FOURIER_LAB_COMPLEX_H

#include <math.h>

/* pi, defined here (not M_PI) so the code builds under strict -std=c11. */
#define FL_PI 3.14159265358979323846

/* A complex number with real and imaginary parts. */
typedef struct {
    double re;
    double im;
} cx_t;

/* Construct a complex number. */
static inline cx_t cx_make(double re, double im) {
    cx_t z;
    z.re = re;
    z.im = im;
    return z;
}

static inline cx_t cx_add(cx_t a, cx_t b) {
    return cx_make(a.re + b.re, a.im + b.im);
}

static inline cx_t cx_sub(cx_t a, cx_t b) {
    return cx_make(a.re - b.re, a.im - b.im);
}

static inline cx_t cx_mul(cx_t a, cx_t b) {
    return cx_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

static inline cx_t cx_scale(cx_t a, double s) {
    return cx_make(a.re * s, a.im * s);
}

static inline cx_t cx_conj(cx_t a) {
    return cx_make(a.re, -a.im);
}

/* Magnitude |z|. */
static inline double cx_abs(cx_t a) {
    return sqrt(a.re * a.re + a.im * a.im);
}

/* Phase arg(z) in (-pi, pi]. */
static inline double cx_arg(cx_t a) {
    return atan2(a.im, a.re);
}

#endif /* FOURIER_LAB_COMPLEX_H */
