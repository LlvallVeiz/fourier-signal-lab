/* signal.h -- Synthesis of test signals.
 *
 * All generators write n samples at sampling rate fs (Hz) into x,
 * starting at t = 0. They are deterministic except for the noise
 * generators, which take an explicit 64-bit seed so every run is
 * reproducible. All functions return 0 on success, -1 on invalid
 * input (NULL pointer, n == 0, fs <= 0, or an out-of-range parameter).
 */

#ifndef FOURIER_LAB_SIGNAL_H
#define FOURIER_LAB_SIGNAL_H

#include <stddef.h>
#include <stdint.h>

/* One sine tone: amp * sin(2*pi*freq*t + phase). */
typedef struct {
    double freq;  /* Hz */
    double amp;
    double phase; /* radians */
} tone_t;

/* Sum of sine tones. tones may be NULL when ntones == 0 (silence). */
int signal_tones(double *x, size_t n, double fs, const tone_t *tones,
                 size_t ntones);

/* Square wave of amplitude amp built from its Fourier series:
 *   (4*amp/pi) * sum_{odd k <= 2*harmonics-1} sin(2*pi*k*freq*t) / k
 * harmonics = number of odd terms (1 gives a pure sine). */
int signal_square(double *x, size_t n, double fs, double freq, double amp,
                  int harmonics);

/* Sawtooth wave of amplitude amp built from its Fourier series:
 *   (2*amp/pi) * sum_{k=1..harmonics} (-1)^(k+1) sin(2*pi*k*freq*t) / k */
int signal_sawtooth(double *x, size_t n, double fs, double freq, double amp,
                    int harmonics);

/* Synthetic ECG: one heartbeat per 60/bpm seconds, each beat a sum of
 * Gaussian bumps (P wave, QRS complex, T wave) in the style of the
 * textbook synthetic ECG used to test QRS detectors such as
 * Pan-Tompkins (1985). The first beat is centred half a period into
 * the record so no beat is clipped at the edge. Optional white
 * Gaussian noise of standard deviation noise_std is added
 * (noise_std == 0 disables it). R-peak amplitude is 1.0 mV. */
int signal_ecg(double *x, size_t n, double fs, double bpm,
               double noise_std, uint64_t seed);

/* White Gaussian noise, mean 0, standard deviation std. */
int signal_noise(double *x, size_t n, double std, uint64_t seed);

/* Linear chirp sweeping from f0 to f1 Hz over the record, amplitude amp. */
int signal_chirp(double *x, size_t n, double fs, double f0, double f1,
                 double amp);

/* Deterministic PRNG (xorshift64*) exposed for tests. */
uint64_t prng_next(uint64_t *state);
double prng_uniform(uint64_t *state); /* uniform in [0, 1) */
double prng_gauss(uint64_t *state);   /* standard normal (Box-Muller) */

#endif /* FOURIER_LAB_SIGNAL_H */
