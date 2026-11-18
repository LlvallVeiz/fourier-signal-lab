/* test_signal.c -- Tests for the signal synthesizers.
 *
 * Hand-computed values used here:
 * - Sine 25 Hz at fs=100 Hz: samples 0,1,2,3 are 0, +1, 0, -1
 *   (quarter-period steps of pi/2).
 * - Square wave, 1 odd harmonic: (4/pi)*sin(2*pi*f*t); at a quarter
 *   period that is 4/pi = 1.2732. With 25 odd harmonics the value at
 *   the quarter period is (4/pi) times the Leibniz partial sum
 *   1 - 1/3 + 1/5 - ... + 1/49 ~= 0.7952 (odd partial sums sit above
 *   pi/4 ~= 0.7854 by ~0.0098), i.e. ~= 1.0127, and the Gibbs
 *   overshoot near the jump is ~= 1.179 (9% of the 2-unit jump).
 * - Sawtooth at a quarter period: (2/pi) times the same Leibniz
 *   partial sum ~= 0.5.
 * - Synthetic ECG at 72 BPM: R spikes of 1.0 mV at t = 0.4167 +
 *   k*0.8333 s; S dips to about -0.25 mV.
 */

#include "test_framework.h"

#include <math.h>
#include <string.h>

#include "complex.h"
#include "signal.h"

void test_signal(void) {
    /* --- Single sine: quarter-period samples are 0, 1, 0, -1. --- */
    {
        double x[8];
        tone_t t = {25.0, 1.0, 0.0};
        CHECK(signal_tones(x, 8, 100.0, &t, 1) == 0);
        CHECK_NEAR(x[0], 0.0, 1e-12);
        CHECK_NEAR(x[1], 1.0, 1e-12);
        CHECK_NEAR(x[2], 0.0, 1e-12);
        CHECK_NEAR(x[3], -1.0, 1e-12);
        CHECK_NEAR(x[4], 0.0, 1e-12);
    }

    /* --- Two tones equal the direct sum, sample by sample. --- */
    {
        double x[64];
        tone_t tones[2] = {{50.0, 1.0, 0.0}, {120.0, 0.5, 0.3}};
        CHECK(signal_tones(x, 64, 1000.0, tones, 2) == 0);
        for (int i = 0; i < 64; i++) {
            const double t = (double)i / 1000.0;
            const double expect = 1.0 * sin(2.0 * FL_PI * 50.0 * t) +
                                  0.5 * sin(2.0 * FL_PI * 120.0 * t + 0.3);
            CHECK_NEAR(x[i], expect, 1e-12);
        }
        /* Zero tones = silence. */
        CHECK(signal_tones(x, 64, 1000.0, NULL, 0) == 0);
        for (int i = 0; i < 64; i++) {
            CHECK_NEAR(x[i], 0.0, 1e-15);
        }
    }

    /* --- Square wave: 1 harmonic is (4/pi)*sine; 25 harmonics show
     * the Gibbs overshoot (~1.179) and ~0.988 at the quarter period. --- */
    {
        double x[1000];
        CHECK(signal_square(x, 1000, 1000.0, 10.0, 1.0, 1) == 0);
        CHECK_NEAR(x[25], 4.0 / FL_PI, 1e-12); /* t = T/4 */
        CHECK_NEAR(x[0], 0.0, 1e-12);

        CHECK(signal_square(x, 1000, 1000.0, 10.0, 1.0, 25) == 0);
        CHECK_NEAR(x[25], 1.0127, 0.005);
        double maxv = 0.0;
        for (int i = 0; i < 500; i++) { /* first half-period */
            if (x[i] > maxv) {
                maxv = x[i];
            }
        }
        CHECK(maxv > 1.15 && maxv < 1.20); /* Gibbs overshoot */
        /* Odd symmetry: x(t + T/2) = -x(t). */
        for (int i = 0; i < 50; i++) {
            CHECK_NEAR(x[i + 50], -x[i], 1e-9);
        }
    }

    /* --- Sawtooth: ~0.5 at a quarter period, 0 at t = 0. --- */
    {
        double x[1000];
        CHECK(signal_sawtooth(x, 1000, 1000.0, 10.0, 1.0, 50) == 0);
        CHECK_NEAR(x[0], 0.0, 1e-12);
        CHECK_NEAR(x[25], 0.5, 0.02);
        /* Rises (roughly) linearly across the first half-period. */
        CHECK(x[10] > 0.0 && x[10] < x[25]);
        CHECK(x[40] > x[25]);
    }

    /* --- Synthetic ECG: 12 R spikes of 1.0 mV in 10 s at 72 BPM. --- */
    {
        double x[3600];
        CHECK(signal_ecg(x, 3600, 360.0, 72.0, 0.0, 42) == 0);
        double maxv = -1e30;
        double minv = 1e30;
        for (int i = 0; i < 3600; i++) {
            if (x[i] > maxv) {
                maxv = x[i];
            }
            if (x[i] < minv) {
                minv = x[i];
            }
        }
        CHECK(maxv > 0.95 && maxv < 1.05);
        CHECK(minv > -0.30 && minv < -0.20); /* S wave */
        /* First R spike at t = 0.5 * (60/72) = 0.4167 s -> sample 150,
         * and it attains (within rounding) the record maximum. */
        CHECK_NEAR(x[150], 1.0, 0.01);
        CHECK(maxv - x[150] < 1e-9);
        /* Count samples above 0.9 mV: ~3 per beat * 12 beats. */
        int hot = 0;
        for (int i = 0; i < 3600; i++) {
            if (x[i] > 0.9) {
                hot++;
            }
        }
        CHECK(hot >= 30 && hot <= 60);
        /* Deterministic: same seed, same signal. */
        double y[3600];
        CHECK(signal_ecg(y, 3600, 360.0, 72.0, 0.0, 42) == 0);
        CHECK(memcmp(x, y, sizeof(x)) == 0);
        /* A different seed with noise differs; noise raises the peak. */
        CHECK(signal_ecg(y, 3600, 360.0, 72.0, 0.05, 7) == 0);
        CHECK(memcmp(x, y, sizeof(x)) != 0);
    }

    /* --- Noise: mean ~0, std ~1, reproducible per seed. --- */
    {
        const size_t n = 100000;
        double x[100000];
        CHECK(signal_noise(x, n, 1.0, 1234) == 0);
        double mean = 0.0;
        double var = 0.0;
        for (size_t i = 0; i < n; i++) {
            mean += x[i];
        }
        mean /= (double)n;
        for (size_t i = 0; i < n; i++) {
            var += (x[i] - mean) * (x[i] - mean);
        }
        var /= (double)n;
        CHECK(fabs(mean) < 0.02);
        CHECK(fabs(sqrt(var) - 1.0) < 0.02);
        double y[100000];
        CHECK(signal_noise(y, n, 1.0, 1234) == 0);
        CHECK(memcmp(x, y, sizeof(x)) == 0);
        /* std = 2 doubles the spread. */
        CHECK(signal_noise(y, n, 2.0, 1234) == 0);
        double mean2 = 0.0;
        for (size_t i = 0; i < n; i++) {
            mean2 += y[i];
        }
        mean2 /= (double)n;
        CHECK(fabs(mean2) < 0.05);
        CHECK(fabs(y[0]) > 0.0); /* non-degenerate */
    }

    /* --- Chirp: starts at 0, stays bounded, rises in frequency. --- */
    {
        double x[1000];
        CHECK(signal_chirp(x, 1000, 1000.0, 10.0, 100.0, 1.0) == 0);
        CHECK_NEAR(x[0], 0.0, 1e-12);
        for (int i = 0; i < 1000; i++) {
            CHECK(fabs(x[i]) <= 1.0 + 1e-12);
        }
        /* Zero crossings get closer together as frequency rises:
         * count crossings in the first vs the last 100 ms. */
        int early = 0, late = 0;
        for (int i = 1; i < 100; i++) {
            if ((x[i - 1] < 0) != (x[i] < 0)) early++;
        }
        for (int i = 901; i < 1000; i++) {
            if ((x[i - 1] < 0) != (x[i] < 0)) late++;
        }
        CHECK(late > early);
    }

    /* --- Invalid inputs are rejected. --- */
    {
        double x[8];
        tone_t t = {10.0, 1.0, 0.0};
        CHECK(signal_tones(NULL, 8, 100.0, &t, 1) == -1);
        CHECK(signal_tones(x, 0, 100.0, &t, 1) == -1);
        CHECK(signal_tones(x, 8, 0.0, &t, 1) == -1);
        CHECK(signal_tones(x, 8, 100.0, NULL, 1) == -1);
        CHECK(signal_square(x, 8, 100.0, 10.0, 1.0, 0) == -1);
        CHECK(signal_sawtooth(x, 8, 100.0, -5.0, 1.0, 5) == -1);
        CHECK(signal_ecg(x, 8, 100.0, 0.0, 0.0, 1) == -1);
        CHECK(signal_ecg(x, 8, 100.0, 72.0, -1.0, 1) == -1);
        CHECK(signal_noise(x, 8, -1.0, 1) == -1);
        CHECK(signal_chirp(x, 8, 100.0, -1.0, 10.0, 1.0) == -1);
    }
}
