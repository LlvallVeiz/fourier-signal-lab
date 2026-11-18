/* signal.c -- Synthesis of test signals (tones, Fourier-series waves,
 * synthetic ECG, noise, chirp).
 *
 * The square and sawtooth waves are built from their Fourier series
 * rather than from a sign() switch, so the Gibbs overshoot at the
 * discontinuities is visible and the harmonic amplitudes follow the
 * textbook 1/k (odd k only for the square wave) law exactly.
 *
 * The ECG model places one Gaussian bump per wave of the heartbeat:
 * P (0.15 mV), Q (-0.10 mV), R (1.00 mV), S (-0.25 mV), T (0.30 mV),
 * with widths chosen so the QRS complex is narrow (~25 ms) and the P
 * and T waves are broad. It is a teaching model, not a physiological
 * simulation: it reproduces the morphology a QRS detector keys on
 * (a sharp R spike on a slow baseline) and nothing more.
 */

#include "signal.h"

#include <math.h>

#include "complex.h"

uint64_t prng_next(uint64_t *state) {
    uint64_t x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * 0x2545F4914F6CDD1DULL;
}

double prng_uniform(uint64_t *state) {
    /* 53 random bits -> double in [0, 1). */
    return (double)(prng_next(state) >> 11) * (1.0 / 9007199254740992.0);
}

double prng_gauss(uint64_t *state) {
    /* Box-Muller: two uniforms -> one standard normal. The second
     * normal is discarded to keep the call stateless. */
    double u1 = prng_uniform(state);
    if (u1 <= 0.0) {
        u1 = 1e-300; /* log(0) guard; unreachable for a 53-bit draw */
    }
    const double u2 = prng_uniform(state);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * FL_PI * u2);
}

int signal_tones(double *x, size_t n, double fs, const tone_t *tones,
                 size_t ntones) {
    if (x == NULL || n == 0 || fs <= 0.0) {
        return -1;
    }
    if (ntones > 0 && tones == NULL) {
        return -1;
    }
    for (size_t i = 0; i < n; i++) {
        const double t = (double)i / fs;
        double v = 0.0;
        for (size_t k = 0; k < ntones; k++) {
            v += tones[k].amp *
                 sin(2.0 * FL_PI * tones[k].freq * t + tones[k].phase);
        }
        x[i] = v;
    }
    return 0;
}

int signal_square(double *x, size_t n, double fs, double freq, double amp,
                  int harmonics) {
    if (x == NULL || n == 0 || fs <= 0.0 || freq <= 0.0 || harmonics < 1) {
        return -1;
    }
    const double scale = 4.0 * amp / FL_PI;
    for (size_t i = 0; i < n; i++) {
        const double t = (double)i / fs;
        double v = 0.0;
        for (int h = 0; h < harmonics; h++) {
            const double k = 2.0 * (double)h + 1.0; /* odd harmonics */
            v += sin(2.0 * FL_PI * k * freq * t) / k;
        }
        x[i] = scale * v;
    }
    return 0;
}

int signal_sawtooth(double *x, size_t n, double fs, double freq, double amp,
                    int harmonics) {
    if (x == NULL || n == 0 || fs <= 0.0 || freq <= 0.0 || harmonics < 1) {
        return -1;
    }
    const double scale = 2.0 * amp / FL_PI;
    for (size_t i = 0; i < n; i++) {
        const double t = (double)i / fs;
        double v = 0.0;
        for (int k = 1; k <= harmonics; k++) {
            const double sign = (k % 2 == 1) ? 1.0 : -1.0; /* (-1)^(k+1) */
            v += sign * sin(2.0 * FL_PI * (double)k * freq * t) / (double)k;
        }
        x[i] = scale * v;
    }
    return 0;
}

/* One Gaussian bump: amp * exp(-0.5 * ((t - centre) / sigma)^2). */
static double gauss_bump(double t, double centre, double sigma, double amp) {
    const double d = (t - centre) / sigma;
    return amp * exp(-0.5 * d * d);
}

int signal_ecg(double *x, size_t n, double fs, double bpm,
               double noise_std, uint64_t seed) {
    if (x == NULL || n == 0 || fs <= 0.0 || bpm <= 0.0 || bpm > 300.0 ||
        noise_std < 0.0) {
        return -1;
    }
    const double period = 60.0 / bpm; /* seconds per beat */
    uint64_t state = seed == 0 ? 0x9E3779B97F4A7C15ULL : seed;
    for (size_t i = 0; i < n; i++) {
        const double t = (double)i / fs;
        /* Time since the nearest beat centre; beats sit at
         * (k + 0.5) * period so the first beat is half a period in. */
        double tb = fmod(t + 0.5 * period, period);
        if (tb > 0.5 * period) {
            tb -= period;
        }
        double v = 0.0;
        v += gauss_bump(tb, -0.16, 0.025, 0.15);  /* P wave  */
        v += gauss_bump(tb, -0.03, 0.010, -0.10); /* Q wave  */
        v += gauss_bump(tb, 0.00, 0.012, 1.00);   /* R spike  */
        v += gauss_bump(tb, 0.03, 0.010, -0.25);  /* S wave  */
        v += gauss_bump(tb, 0.28, 0.050, 0.30);   /* T wave  */
        if (noise_std > 0.0) {
            v += noise_std * prng_gauss(&state);
        }
        x[i] = v;
    }
    return 0;
}

int signal_noise(double *x, size_t n, double std, uint64_t seed) {
    if (x == NULL || n == 0 || std < 0.0) {
        return -1;
    }
    uint64_t state = seed == 0 ? 0x9E3779B97F4A7C15ULL : seed;
    for (size_t i = 0; i < n; i++) {
        x[i] = std * prng_gauss(&state);
    }
    return 0;
}

int signal_chirp(double *x, size_t n, double fs, double f0, double f1,
                 double amp) {
    if (x == NULL || n == 0 || fs <= 0.0 || f0 < 0.0 || f1 < 0.0) {
        return -1;
    }
    const double duration = (double)n / fs;
    for (size_t i = 0; i < n; i++) {
        const double t = (double)i / fs;
        /* Phase integrates the linearly rising frequency:
         * phi(t) = 2*pi*(f0*t + (f1-f0)*t^2 / (2*T)). */
        const double phase =
            2.0 * FL_PI * (f0 * t + (f1 - f0) * t * t / (2.0 * duration));
        x[i] = amp * sin(phase);
    }
    return 0;
}
