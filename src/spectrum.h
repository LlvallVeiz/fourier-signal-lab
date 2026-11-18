/* spectrum.h -- Reading a spectrum: single-sided amplitude/phase,
 * peak detection with parabolic frequency refinement, Parseval checks,
 * spectral-leakage measurement, and a simplified Pan-Tompkins QRS
 * detector for the synthetic ECG.
 */

#ifndef FOURIER_LAB_SPECTRUM_H
#define FOURIER_LAB_SPECTRUM_H

#include <stddef.h>

#include "complex.h"
#include "window.h"

/* One detected spectral peak. */
typedef struct {
    size_t bin;       /* FFT bin index */
    double freq;      /* Hz, parabolically refined */
    double magnitude; /* amplitude (same units as the signal) */
    double phase;     /* radians */
} peak_t;

/* Single-sided amplitude spectrum of a real signal.
 * X is the N-point FFT of the windowed signal; mag has N/2+1 entries
 * (bins 0..N/2). Amplitudes are corrected for the window's coherent
 * gain, so a tone of amplitude A reads A at its bin: mag[k] =
 * (2/N) * |X[k]| * correction, with the factor 2 dropped for the DC
 * and Nyquist bins (which have no mirror image). */
void spectrum_magnitude(const cx_t *X, size_t n, window_t w, double *mag);

/* Phase spectrum (radians), same layout as spectrum_magnitude. */
void spectrum_phase(const cx_t *X, size_t n, double *phase);

/* Detect peaks in a magnitude spectrum of nbins bins with bin width
 * df Hz. A peak is a local maximum above min_height; accepted peaks
 * are at least min_distance bins apart (strongest first). Results are
 * sorted by frequency ascending. Returns the number of peaks written
 * (<= max_peaks), or -1 on invalid input. */
int spectrum_find_peaks(const double *mag, size_t nbins, double df,
                        double min_height, size_t min_distance,
                        peak_t *peaks, size_t max_peaks);

/* Parabolic interpolation of a peak at bin k (1 <= k <= nbins-2):
 * returns the fractional bin offset in (-1, 1). The parabola is
 * fitted to the *log* magnitude, which is the standard refinement
 * (Smith, "Spectral Audio Signal Processing"): it is near-exact for
 * Hann-windowed tones and much better than plain-magnitude fitting
 * for other windows. Falls back to plain magnitudes when any of the
 * three is zero. */
double spectrum_parabolic_delta(const double *mag, size_t k);

/* Refined peak frequency in Hz for a peak at bin k. */
double spectrum_parabolic_freq(const double *mag, size_t k, double df);

/* Parseval: sum(x^2) and (1/N) * sum(|X|^2) must agree for a correct
 * transform pair. */
double spectrum_parseval_time(const double *x, size_t n);
double spectrum_parseval_freq(const cx_t *X, size_t n);

/* Fraction of spectral energy (|X|^2) lying OUTSIDE +-spread bins of
 * the dominant peak: 0 for a perfectly concentrated tone, near 1 for
 * a tone smeared across the spectrum by leakage. */
double spectrum_leakage(const cx_t *X, size_t n, size_t spread);

/* Simplified Pan-Tompkins QRS detection on a raw ECG (Pan & Tompkins,
 * 1985): 5-point derivative -> squaring -> 150 ms moving-window
 * integration -> threshold at 25% of the maximum, with a 200 ms
 * refractory period, and a refinement step that takes the largest raw
 * sample within +-100 ms of each hit (the R spike). Writes sample
 * indices of detected R peaks. Returns the number found
 * (<= max_peaks), or -1 on invalid input. This is the textbook
 * pipeline reduced to its essentials; it is validated here only on
 * the synthetic ECG from signal_ecg(). */
int spectrum_qrs_detect(const double *x, size_t n, double fs,
                        size_t *r_peaks, size_t max_peaks);

/* Heart rate in BPM from detected R-peak indices (mean RR interval).
 * Returns -1.0 when fewer than two peaks are given. */
double spectrum_bpm(const size_t *r_peaks, size_t count, double fs);

#endif /* FOURIER_LAB_SPECTRUM_H */
