/* test_spectrum.c -- Tests for spectrum reading, peak detection,
 * parabolic interpolation, Parseval, leakage, and QRS detection.
 *
 * Setup used throughout: fs = 1024 Hz, N = 1024, so the bin width is
 * exactly 1 Hz and tones at integer frequencies sit exactly on bins.
 */

#include "test_framework.h"

#include <math.h>
#include <stdlib.h>

#include "complex.h"
#include "fft.h"
#include "signal.h"
#include "spectrum.h"
#include "window.h"

/* FFT the real signal x (length n, already a power of two) with the
 * given window, into spec (length n). */
static void fft_windowed(const double *x, size_t n, window_t w, cx_t *spec) {
    double *win = (double *)malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        win[i] = x[i];
    }
    window_apply(w, win, n);
    for (size_t i = 0; i < n; i++) {
        spec[i] = cx_make(win[i], 0.0);
    }
    free(win);
    fft_forward(spec, n);
}

void test_spectrum(void) {
    const size_t n = 1024;
    const double fs = 1024.0;
    const size_t nbins = n / 2 + 1;

    /* --- Magnitude of an on-bin tone reads its amplitude, with any
     * window, once the amplitude correction is applied. --- */
    {
        double x[1024];
        tone_t t = {50.0, 2.0, 0.0};
        CHECK(signal_tones(x, n, fs, &t, 1) == 0);
        cx_t spec[1024];
        double mag[513];
        fft_windowed(x, n, WINDOW_RECTANGULAR, spec);
        spectrum_magnitude(spec, n, WINDOW_RECTANGULAR, mag);
        CHECK_NEAR(mag[50], 2.0, 1e-9);
        CHECK_NEAR(mag[0], 0.0, 1e-9);
        CHECK_NEAR(mag[120], 0.0, 1e-9);
        fft_windowed(x, n, WINDOW_HANN, spec);
        spectrum_magnitude(spec, n, WINDOW_HANN, mag);
        CHECK_NEAR(mag[50], 2.0, 1e-9);
        fft_windowed(x, n, WINDOW_HAMMING, spec);
        spectrum_magnitude(spec, n, WINDOW_HAMMING, mag);
        CHECK_NEAR(mag[50], 2.0, 1e-9);
        fft_windowed(x, n, WINDOW_BLACKMAN, spec);
        spectrum_magnitude(spec, n, WINDOW_BLACKMAN, mag);
        CHECK_NEAR(mag[50], 2.0, 1e-9);
    }

    /* --- Phase of a sine at its bin is -pi/2 (sine = -i * ... ). --- */
    {
        double x[1024];
        tone_t t = {50.0, 1.0, 0.0};
        CHECK(signal_tones(x, n, fs, &t, 1) == 0);
        cx_t spec[1024];
        double phase[513];
        fft_windowed(x, n, WINDOW_RECTANGULAR, spec);
        spectrum_phase(spec, n, phase);
        CHECK_NEAR(phase[50], -FL_PI / 2.0, 1e-9);
        /* A cosine of the same frequency has phase 0. */
        tone_t c = {50.0, 1.0, FL_PI / 2.0};
        CHECK(signal_tones(x, n, fs, &c, 1) == 0);
        fft_windowed(x, n, WINDOW_RECTANGULAR, spec);
        spectrum_phase(spec, n, phase);
        CHECK_NEAR(phase[50], 0.0, 1e-9);
    }

    /* --- Peak detection finds both tones of a two-tone signal. --- */
    {
        double x[1024];
        tone_t tones[2] = {{50.0, 1.0, 0.0}, {120.0, 0.5, 0.0}};
        CHECK(signal_tones(x, n, fs, tones, 2) == 0);
        cx_t spec[1024];
        double mag[513];
        fft_windowed(x, n, WINDOW_HANN, spec);
        spectrum_magnitude(spec, n, WINDOW_HANN, mag);
        peak_t peaks[5];
        const int found =
            spectrum_find_peaks(mag, nbins, 1.0, 0.05, 2, peaks, 5);
        CHECK(found == 2);
        CHECK_NEAR(peaks[0].freq, 50.0, 0.01);
        CHECK_NEAR(peaks[0].magnitude, 1.0, 0.02);
        CHECK_NEAR(peaks[1].freq, 120.0, 0.01);
        CHECK_NEAR(peaks[1].magnitude, 0.5, 0.02);
        /* A high threshold hides the weaker tone. */
        const int one =
            spectrum_find_peaks(mag, nbins, 1.0, 0.8, 2, peaks, 5);
        CHECK(one == 1);
        CHECK_NEAR(peaks[0].freq, 50.0, 0.01);
    }

    /* --- Parabolic interpolation refines an off-bin tone. --- */
    {
        double x[1024];
        tone_t t = {50.3, 1.0, 0.0};
        CHECK(signal_tones(x, n, fs, &t, 1) == 0);
        cx_t spec[1024];
        double mag[513];
        fft_windowed(x, n, WINDOW_RECTANGULAR, spec);
        spectrum_magnitude(spec, n, WINDOW_RECTANGULAR, mag);
        peak_t peaks[5];
        const int found =
            spectrum_find_peaks(mag, nbins, 1.0, 0.05, 2, peaks, 5);
        CHECK(found == 1);
        /* Nearest bin is 50; interpolation must pull toward 50.3.
         * Log-magnitude interpolation is crude for a rectangular
         * window (measured error 0.17 bins), so the tolerance is
         * loose here; the Hann case below is the precise one. */
        CHECK(peaks[0].bin == 50);
        CHECK(peaks[0].freq > 50.0 && peaks[0].freq < 50.6);
        CHECK_NEAR(peaks[0].freq, 50.3, 0.25);
        /* With a Hann window the interpolation is tight (measured
         * error 0.016 bins). */
        fft_windowed(x, n, WINDOW_HANN, spec);
        spectrum_magnitude(spec, n, WINDOW_HANN, mag);
        const int found_h =
            spectrum_find_peaks(mag, nbins, 1.0, 0.05, 2, peaks, 5);
        CHECK(found_h == 1);
        CHECK_NEAR(peaks[0].freq, 50.3, 0.05);
    }

    /* --- Parseval helpers agree for a random signal. --- */
    {
        double x[1024];
        unsigned long state = 31337;
        for (size_t i = 0; i < n; i++) {
            state = state * 6364136223846793005UL + 1442695040888963407UL;
            x[i] = (double)((state >> 11) & 0x1FFFFFFFFFFFFFULL) /
                       4503599627370496.0 * 2.0 -
                   1.0;
        }
        cx_t spec[1024];
        fft_windowed(x, n, WINDOW_RECTANGULAR, spec);
        const double lhs = spectrum_parseval_time(x, n);
        const double rhs = spectrum_parseval_freq(spec, n);
        CHECK(fabs(lhs - rhs) / lhs < 1e-12);
    }

    /* --- Leakage: an off-bin tone smears with rectangular, stays
     * concentrated with Hann. --- */
    {
        double x[1024];
        tone_t t = {55.3, 1.0, 0.0};
        CHECK(signal_tones(x, n, fs, &t, 1) == 0);
        cx_t spec[1024];
        fft_windowed(x, n, WINDOW_RECTANGULAR, spec);
        const double leak_rect = spectrum_leakage(spec, n, 2);
        fft_windowed(x, n, WINDOW_HANN, spec);
        const double leak_hann = spectrum_leakage(spec, n, 2);
        /* Measured for a 55.3 Hz tone (N=1024): rect 0.0531,
         * Hann 0.0002 -- a factor of ~200 less leakage energy. */
        CHECK(leak_rect > 0.04);
        CHECK(leak_hann < 0.01);
        CHECK(leak_hann < leak_rect / 10.0);
        /* An on-bin tone does not leak at all, either way. */
        tone_t on = {55.0, 1.0, 0.0};
        CHECK(signal_tones(x, n, fs, &on, 1) == 0);
        fft_windowed(x, n, WINDOW_RECTANGULAR, spec);
        CHECK(spectrum_leakage(spec, n, 2) < 1e-9);
    }

    /* --- Nyquist / aliasing: a 700 Hz tone sampled at 1024 Hz
     * appears at 1024 - 700 = 324 Hz. --- */
    {
        double x[1024];
        tone_t t = {700.0, 1.0, 0.0};
        CHECK(signal_tones(x, n, fs, &t, 1) == 0);
        cx_t spec[1024];
        double mag[513];
        fft_windowed(x, n, WINDOW_RECTANGULAR, spec);
        spectrum_magnitude(spec, n, WINDOW_RECTANGULAR, mag);
        peak_t peaks[5];
        const int found =
            spectrum_find_peaks(mag, nbins, 1.0, 0.05, 2, peaks, 5);
        CHECK(found == 1);
        CHECK_NEAR(peaks[0].freq, 324.0, 0.01);
        CHECK_NEAR(peaks[0].magnitude, 1.0, 0.02);
    }

    /* --- Square-wave harmonics follow the 1/k law (odd k only). --- */
    {
        double x[1024];
        CHECK(signal_square(x, n, fs, 10.0, 1.0, 25) == 0);
        cx_t spec[1024];
        double mag[513];
        fft_windowed(x, n, WINDOW_HANN, spec);
        spectrum_magnitude(spec, n, WINDOW_HANN, mag);
        /* Fundamental 4/pi = 1.2732; 3rd harmonic 1/3 of it; 5th 1/5. */
        CHECK_NEAR(mag[10], 4.0 / FL_PI, 0.02);
        CHECK_NEAR(mag[30], 4.0 / FL_PI / 3.0, 0.02);
        CHECK_NEAR(mag[50], 4.0 / FL_PI / 5.0, 0.02);
        /* Even harmonics are (nearly) absent. */
        CHECK(mag[20] < 0.01);
        CHECK(mag[40] < 0.01);
    }

    /* --- QRS detection on the synthetic ECG: 12 beats, 72 BPM. --- */
    {
        const size_t m = 3600;
        const double ecg_fs = 360.0;
        double x[3600];
        CHECK(signal_ecg(x, m, ecg_fs, 72.0, 0.0, 42) == 0);
        size_t peaks[64];
        const int count = spectrum_qrs_detect(x, m, ecg_fs, peaks, 64);
        CHECK(count == 12);
        if (count == 12) {
            /* True R spikes sit at samples 150 + 300*k (beats at
             * 0.4167 + k*0.8333 s); the refinement step must land
             * within 2 samples of them. */
            for (int i = 0; i < 12; i++) {
                CHECK(peaks[i] >= (size_t)(150 + 300 * i - 2) &&
                      peaks[i] <= (size_t)(150 + 300 * i + 2));
            }
            CHECK_NEAR(spectrum_bpm(peaks, 12, ecg_fs), 72.0, 0.5);
        }
        /* Noisy ECG (std 0.05 mV) is still detected. */
        CHECK(signal_ecg(x, m, ecg_fs, 72.0, 0.05, 7) == 0);
        const int count_noisy = spectrum_qrs_detect(x, m, ecg_fs, peaks, 64);
        CHECK(count_noisy == 12);
        if (count_noisy == 12) {
            CHECK_NEAR(spectrum_bpm(peaks, 12, ecg_fs), 72.0, 0.5);
        }
        /* Silence detects nothing; too few peaks give no BPM. */
        for (size_t i = 0; i < m; i++) {
            x[i] = 0.0;
        }
        CHECK(spectrum_qrs_detect(x, m, ecg_fs, peaks, 64) == 0);
        CHECK_NEAR(spectrum_bpm(peaks, 1, ecg_fs), -1.0, 1e-9);
        CHECK(spectrum_qrs_detect(NULL, m, ecg_fs, peaks, 64) == -1);
        CHECK(spectrum_qrs_detect(x, 8, ecg_fs, peaks, 64) == -1);
    }

    /* --- Invalid inputs are rejected. --- */
    {
        double mag[513];
        peak_t peaks[5];
        CHECK(spectrum_find_peaks(NULL, nbins, 1.0, 0.05, 2, peaks, 5) == -1);
        CHECK(spectrum_find_peaks(mag, 2, 1.0, 0.05, 2, peaks, 5) == -1);
        CHECK(spectrum_find_peaks(mag, nbins, 1.0, 0.05, 2, peaks, 0) == -1);
    }
}
