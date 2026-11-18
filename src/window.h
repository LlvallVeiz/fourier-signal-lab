/* window.h -- Data windows for spectral analysis.
 *
 * A finite record is implicitly periodic to the DFT, so a tone that
 * does not complete an integer number of cycles in the record leaks
 * energy into neighbouring bins. Multiplying the record by a window
 * tapers its ends to zero and trades main-lobe width against
 * sidelobe level (Harris, 1978). All windows here are the *periodic*
 * (DFT) form w[n], n = 0..N-1, so that a tone exactly on a bin is
 * recovered exactly after amplitude correction.
 *
 * Figures of merit (coherent gain, ENBW, amplitude correction) are
 * computed numerically from the coefficients, so they are exact for
 * the N actually used; the highest-sidelobe figures are the textbook
 * values tabulated by Harris (1978, Table I).
 */

#ifndef FOURIER_LAB_WINDOW_H
#define FOURIER_LAB_WINDOW_H

#include <stddef.h>

typedef enum {
    WINDOW_RECTANGULAR = 0,
    WINDOW_HANN,
    WINDOW_HAMMING,
    WINDOW_BLACKMAN,
    WINDOW_COUNT
} window_t;

/* Parse a window name ("rect", "rectangular", "hann", "hanning",
 * "hamming", "blackman"; case-insensitive). Returns 0 on success. */
int window_parse(const char *name, window_t *out);

/* Short display name, e.g. "hann". Never NULL. */
const char *window_name(window_t w);

/* Coefficient w[n] for 0 <= n < N. Returns 0.0 for invalid input. */
double window_coef(window_t w, size_t n, size_t N);

/* Multiply x[0..n-1] by the window in place. */
void window_apply(window_t w, double *x, size_t n);

/* Coherent gain = (1/N) * sum(w): the factor by which a coherent
 * (on-bin) tone's amplitude is scaled by the window. */
double window_coherent_gain(window_t w, size_t n);

/* Equivalent noise bandwidth in bins: N * sum(w^2) / sum(w)^2.
 * Rectangular = 1.0 by definition. */
double window_enbw(window_t w, size_t n);

/* Amplitude correction = 1 / coherent gain. Multiply a windowed
 * spectrum by this to recover true tone amplitudes. */
double window_amplitude_correction(window_t w, size_t n);

/* Highest sidelobe level in dB (negative), Harris (1978) Table I:
 * rectangular -13, Hann -31, Hamming -43, Blackman -58. */
double window_highest_sidelobe_db(window_t w);

#endif /* FOURIER_LAB_WINDOW_H */
