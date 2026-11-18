/* main.c -- fourier-lab command line interface.
 *
 *   fourier-lab fft         --input in.csv --output spec.csv [--n N]
 *                            [--window hann] [--fs Hz]
 *   fourier-lab analyze     --input in.csv [--fs Hz] [--window hann]
 *                            [--n N] [--min-height H] [--top K]
 *   fourier-lab synthesize  --type tones|square|sawtooth|ecg|noise|chirp
 *                            --output out.csv --fs Hz --duration S [...]
 *   fourier-lab reconstruct --input in.csv [--fs Hz]
 *   fourier-lab qrs         --input in.csv [--fs Hz]
 *
 * Signals are plain-text CSV files: one sample per line, with an
 * optional "# fs=<rate>" header. See README.md for worked examples.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "complex.h"
#include "fft.h"
#include "io.h"
#include "signal.h"
#include "spectrum.h"
#include "window.h"

#define FOURIER_LAB_VERSION "1.0.0"

static void print_usage(FILE *out) {
    fprintf(out,
        "fourier-lab " FOURIER_LAB_VERSION
        " -- Fourier analysis lab (C, no dependencies)\n"
        "\n"
        "Usage:\n"
        "  fourier-lab fft --input IN.csv --output OUT.csv [--n N]\n"
        "                  [--window hann] [--fs Hz]\n"
        "  fourier-lab analyze --input IN.csv [--fs Hz] [--window hann]\n"
        "                  [--n N] [--min-height H] [--top K]\n"
        "  fourier-lab synthesize --type T --output OUT.csv --fs Hz\n"
        "                  --duration S [type options]\n"
        "  fourier-lab reconstruct --input IN.csv [--fs Hz]\n"
        "  fourier-lab qrs --input IN.csv [--fs Hz]\n"
        "  fourier-lab --help | --version\n"
        "\n"
        "Windows: rectangular, hann, hamming, blackman (default: hann).\n"
        "Synthesize types: tones, square, sawtooth, ecg, noise, chirp.\n"
        "  tones:    --tones \"50:1.0,120:0.5\"  (freq:amp pairs)\n"
        "  square:   --freq Hz --amp A --harmonics K (odd terms)\n"
        "  sawtooth: --freq Hz --amp A --harmonics K\n"
        "  ecg:      --bpm B [--noise-std S] [--seed N]\n"
        "  noise:    --noise-std S [--seed N]\n"
        "  chirp:    --f0 Hz --f1 Hz --amp A\n"
        "\n"
        "Signal CSV: one sample per line, optional \"# fs=<rate>\" header.\n");
}

/* ---- small option parser ------------------------------------------- */

/* Find the value of "--name value" or "--name=value" in argv.
 * Returns NULL when the option is absent. */
static const char *opt_value(int argc, char **argv, const char *name) {
    const size_t len = strlen(name);
    for (int i = 0; i < argc; i++) {
        if (strncmp(argv[i], "--", 2) != 0) {
            continue;
        }
        const char *key = argv[i] + 2;
        if (strncmp(key, name, len) == 0 && key[len] == '=') {
            return key + len + 1;
        }
        if (strcmp(key, name) == 0 && i + 1 < argc) {
            return argv[i + 1];
        }
    }
    return NULL;
}

static int opt_double(int argc, char **argv, const char *name, double *out) {
    const char *v = opt_value(argc, argv, name);
    if (v == NULL) {
        return 0; /* absent: caller keeps its default */
    }
    char *end = NULL;
    const double d = strtod(v, &end);
    if (end == v || *end != '\0') {
        return -1;
    }
    *out = d;
    return 1;
}

static int opt_long(int argc, char **argv, const char *name, long *out) {
    const char *v = opt_value(argc, argv, name);
    if (v == NULL) {
        return 0;
    }
    char *end = NULL;
    const long d = strtol(v, &end, 10);
    if (end == v || *end != '\0') {
        return -1;
    }
    *out = d;
    return 1;
}

/* Load a signal, resolving fs from --fs, else the file header. */
static int load_signal(int argc, char **argv, const char *path,
                       double **x, size_t *n, double *fs) {
    if (io_read_signal_csv(path, x, n, fs) != 0) {
        fprintf(stderr, "fourier-lab: cannot read signal file '%s'\n", path);
        return -1;
    }
    double fs_opt = 0.0;
    if (opt_double(argc, argv, "fs", &fs_opt) < 0) {
        fprintf(stderr, "fourier-lab: invalid --fs value\n");
        free(*x);
        return -1;
    }
    if (fs_opt > 0.0) {
        *fs = fs_opt;
    }
    if (*fs <= 0.0) {
        fprintf(stderr,
                "fourier-lab: sampling rate unknown -- pass --fs or add "
                "a '# fs=<rate>' header to the file\n");
        free(*x);
        return -1;
    }
    if (*n == 0) {
        fprintf(stderr, "fourier-lab: signal file '%s' is empty\n", path);
        free(*x);
        return -1;
    }
    return 0;
}

/* Zero-pad or truncate a real signal to length n and FFT it. The
 * caller owns the returned spectrum (free with free()). */
static cx_t *fft_of_signal(const double *x, size_t len, size_t n,
                           window_t w) {
    double *win = (double *)malloc(n * sizeof(double));
    cx_t *spec = (cx_t *)malloc(n * sizeof(cx_t));
    if (win == NULL || spec == NULL) {
        free(win);
        free(spec);
        return NULL;
    }
    for (size_t i = 0; i < n; i++) {
        win[i] = i < len ? x[i] : 0.0;
    }
    window_apply(w, win, n);
    for (size_t i = 0; i < n; i++) {
        spec[i] = cx_make(win[i], 0.0);
    }
    free(win);
    if (fft_forward(spec, n) != 0) {
        free(spec);
        return NULL;
    }
    return spec;
}

/* ---- subcommands ---------------------------------------------------- */

static int cmd_fft(int argc, char **argv) {
    const char *input = opt_value(argc, argv, "input");
    const char *output = opt_value(argc, argv, "output");
    if (input == NULL || output == NULL) {
        fprintf(stderr, "fourier-lab fft: --input and --output required\n");
        return 2;
    }
    window_t w = WINDOW_HANN;
    const char *wname = opt_value(argc, argv, "window");
    if (wname != NULL && window_parse(wname, &w) != 0) {
        fprintf(stderr, "fourier-lab fft: unknown window '%s'\n", wname);
        return 2;
    }
    long n_opt = 0;
    if (opt_long(argc, argv, "n", &n_opt) < 0 || n_opt < 0) {
        fprintf(stderr, "fourier-lab fft: invalid --n value\n");
        return 2;
    }
    double *x = NULL;
    size_t len = 0;
    double fs = 0.0;
    if (load_signal(argc, argv, input, &x, &len, &fs) != 0) {
        return 1;
    }
    size_t n = (size_t)n_opt;
    if (n == 0) {
        n = fft_is_pow2(len) ? len : fft_next_pow2(len);
    }
    if (!fft_is_pow2(n)) {
        fprintf(stderr, "fourier-lab fft: N=%zu is not a power of two\n", n);
        free(x);
        return 2;
    }
    cx_t *spec = fft_of_signal(x, len, n, w);
    if (spec == NULL) {
        fprintf(stderr, "fourier-lab fft: transform failed\n");
        free(x);
        return 1;
    }
    const size_t nbins = n / 2 + 1;
    double *freq = (double *)malloc(nbins * sizeof(double));
    double *mag = (double *)malloc(nbins * sizeof(double));
    double *phase = (double *)malloc(nbins * sizeof(double));
    if (freq == NULL || mag == NULL || phase == NULL) {
        fprintf(stderr, "fourier-lab fft: out of memory\n");
        free(x);
        free(spec);
        free(freq);
        free(mag);
        free(phase);
        return 1;
    }
    const double df = fs / (double)n;
    for (size_t k = 0; k < nbins; k++) {
        freq[k] = (double)k * df;
    }
    spectrum_magnitude(spec, n, w, mag);
    spectrum_phase(spec, n, phase);
    if (io_write_spectrum_csv(output, freq, mag, phase, nbins, fs, n,
                              window_name(w)) != 0) {
        fprintf(stderr, "fourier-lab fft: cannot write '%s'\n", output);
        free(x);
        free(spec);
        free(freq);
        free(mag);
        free(phase);
        return 1;
    }
    peak_t peaks[5];
    const int npeaks = spectrum_find_peaks(mag, nbins, df, 0.0, 2, peaks, 5);
    printf("FFT: N=%zu fs=%.17g Hz df=%.6g Hz window=%s\n", n, fs, df,
           window_name(w));
    printf("spectrum written to %s (%zu bins)\n", output, nbins);
    if (npeaks > 0) {
        printf("top peaks:\n");
        for (int i = 0; i < npeaks; i++) {
            printf("  %8.3f Hz  magnitude %.6g\n", peaks[i].freq,
                   peaks[i].magnitude);
        }
    }
    free(x);
    free(spec);
    free(freq);
    free(mag);
    free(phase);
    return 0;
}

static int cmd_analyze(int argc, char **argv) {
    const char *input = opt_value(argc, argv, "input");
    if (input == NULL) {
        fprintf(stderr, "fourier-lab analyze: --input required\n");
        return 2;
    }
    window_t w = WINDOW_HANN;
    const char *wname = opt_value(argc, argv, "window");
    if (wname != NULL && window_parse(wname, &w) != 0) {
        fprintf(stderr, "fourier-lab analyze: unknown window '%s'\n", wname);
        return 2;
    }
    long n_opt = 0;
    if (opt_long(argc, argv, "n", &n_opt) < 0 || n_opt < 0) {
        fprintf(stderr, "fourier-lab analyze: invalid --n value\n");
        return 2;
    }
    double min_height = 0.05;
    if (opt_double(argc, argv, "min-height", &min_height) < 0 ||
        min_height < 0.0) {
        fprintf(stderr, "fourier-lab analyze: invalid --min-height value\n");
        return 2;
    }
    long top = 5;
    if (opt_long(argc, argv, "top", &top) < 0 || top < 1 || top > 64) {
        fprintf(stderr, "fourier-lab analyze: invalid --top value\n");
        return 2;
    }
    double *x = NULL;
    size_t len = 0;
    double fs = 0.0;
    if (load_signal(argc, argv, input, &x, &len, &fs) != 0) {
        return 1;
    }
    size_t n = (size_t)n_opt;
    if (n == 0) {
        n = fft_is_pow2(len) ? len : fft_next_pow2(len);
    }
    if (!fft_is_pow2(n)) {
        fprintf(stderr, "fourier-lab analyze: N=%zu is not a power of two\n",
                n);
        free(x);
        return 2;
    }
    cx_t *spec = fft_of_signal(x, len, n, w);
    if (spec == NULL) {
        fprintf(stderr, "fourier-lab analyze: transform failed\n");
        free(x);
        return 1;
    }
    const size_t nbins = n / 2 + 1;
    const double df = fs / (double)n;
    double *mag = (double *)malloc(nbins * sizeof(double));
    if (mag == NULL) {
        fprintf(stderr, "fourier-lab analyze: out of memory\n");
        free(x);
        free(spec);
        return 1;
    }
    spectrum_magnitude(spec, n, w, mag);

    printf("analyze: %s\n", input);
    printf("  N=%zu  fs=%.17g Hz  df=%.6g Hz  Nyquist=%.6g Hz  window=%s\n",
           n, fs, df, fs / 2.0, window_name(w));
    printf("  window: coherent gain %.4f, ENBW %.4f bins, amplitude "
           "correction %.4f, highest sidelobe %.0f dB\n",
           window_coherent_gain(w, n), window_enbw(w, n),
           window_amplitude_correction(w, n),
           window_highest_sidelobe_db(w));

    peak_t *peaks = (peak_t *)malloc((size_t)top * sizeof(peak_t));
    if (peaks == NULL) {
        fprintf(stderr, "fourier-lab analyze: out of memory\n");
        free(x);
        free(spec);
        free(mag);
        return 1;
    }
    const int npeaks =
        spectrum_find_peaks(mag, nbins, df, min_height, 2, peaks, (size_t)top);
    if (npeaks <= 0) {
        printf("  no peaks above %.6g\n", min_height);
    } else {
        printf("  peaks (freq, magnitude):\n");
        for (int i = 0; i < npeaks; i++) {
            printf("    %10.3f Hz  %.6g\n", peaks[i].freq, peaks[i].magnitude);
        }
    }

    /* Parseval check on the windowed signal (uses the raw FFT). */
    double *win = (double *)malloc(n * sizeof(double));
    if (win != NULL) {
        for (size_t i = 0; i < n; i++) {
            win[i] = i < len ? x[i] : 0.0;
        }
        window_apply(w, win, n);
        const double lhs = spectrum_parseval_time(win, n);
        const double rhs = spectrum_parseval_freq(spec, n);
        const double rel = lhs > 0.0 ? fabs(lhs - rhs) / lhs : 0.0;
        printf("  Parseval: sum(x^2)=%.6g  (1/N)sum(|X|^2)=%.6g  "
               "rel. diff=%.3g\n",
               lhs, rhs, rel);
        free(win);
    }
    const double leak = spectrum_leakage(spec, n, 2);
    printf("  leakage: %.4f of spectral energy lies outside +-2 bins of "
           "the dominant peak\n",
           leak);
    if (npeaks > 0 && peaks[0].freq > 0.45 * fs) {
        printf("  note: dominant peak near Nyquist (%.3f Hz) -- it may be "
               "an alias of a higher frequency\n",
               peaks[0].freq);
    }
    printf("  frequency resolution (bin width): %.6g Hz -- two tones "
           "closer than this cannot be separated\n",
           df);
    free(x);
    free(spec);
    free(mag);
    free(peaks);
    return 0;
}

static int cmd_reconstruct(int argc, char **argv) {
    const char *input = opt_value(argc, argv, "input");
    if (input == NULL) {
        fprintf(stderr, "fourier-lab reconstruct: --input required\n");
        return 2;
    }
    double *x = NULL;
    size_t len = 0;
    double fs = 0.0;
    if (load_signal(argc, argv, input, &x, &len, &fs) != 0) {
        return 1;
    }
    const size_t n = fft_is_pow2(len) ? len : fft_next_pow2(len);
    cx_t *spec = (cx_t *)malloc(n * sizeof(cx_t));
    if (spec == NULL) {
        fprintf(stderr, "fourier-lab reconstruct: out of memory\n");
        free(x);
        return 1;
    }
    for (size_t i = 0; i < n; i++) {
        spec[i] = cx_make(i < len ? x[i] : 0.0, 0.0);
    }
    if (fft_forward(spec, n) != 0 || fft_inverse(spec, n) != 0) {
        fprintf(stderr, "fourier-lab reconstruct: transform failed\n");
        free(x);
        free(spec);
        return 1;
    }
    double max_err = 0.0;
    double sq_err = 0.0;
    double max_val = 0.0;
    for (size_t i = 0; i < len; i++) {
        const double err = fabs(spec[i].re - x[i]);
        if (err > max_err) {
            max_err = err;
        }
        sq_err += err * err;
        if (fabs(x[i]) > max_val) {
            max_val = fabs(x[i]);
        }
    }
    const double rms_err = sqrt(sq_err / (double)len);
    printf("reconstruct: %s (N=%zu, fs=%.17g Hz)\n", input, n, fs);
    printf("  FFT -> IFFT round trip over %zu samples\n", len);
    printf("  max abs error: %.6g\n", max_err);
    printf("  RMS error:     %.6g\n", rms_err);
    printf("  max |x|:       %.6g  (relative max error %.6g)\n", max_val,
           max_val > 0.0 ? max_err / max_val : 0.0);
    free(x);
    free(spec);
    return 0;
}

static int cmd_qrs(int argc, char **argv) {
    const char *input = opt_value(argc, argv, "input");
    if (input == NULL) {
        fprintf(stderr, "fourier-lab qrs: --input required\n");
        return 2;
    }
    double *x = NULL;
    size_t len = 0;
    double fs = 0.0;
    if (load_signal(argc, argv, input, &x, &len, &fs) != 0) {
        return 1;
    }
    size_t *peaks = (size_t *)malloc(len * sizeof(size_t));
    if (peaks == NULL) {
        fprintf(stderr, "fourier-lab qrs: out of memory\n");
        free(x);
        return 1;
    }
    const int count = spectrum_qrs_detect(x, len, fs, peaks, len);
    if (count < 0) {
        fprintf(stderr, "fourier-lab qrs: detection failed\n");
        free(x);
        free(peaks);
        return 1;
    }
    printf("qrs: %s (fs=%.17g Hz, %zu samples)\n", input, fs, len);
    printf("  R peaks detected: %d\n", count);
    if (count > 0) {
        printf("  R times (s):");
        for (int i = 0; i < count; i++) {
            printf(" %.4f", (double)peaks[i] / fs);
        }
        printf("\n");
    }
    if (count >= 2) {
        printf("  heart rate: %.2f BPM\n", spectrum_bpm(peaks, (size_t)count, fs));
    } else {
        printf("  heart rate: n/a (need at least 2 beats)\n");
    }
    free(x);
    free(peaks);
    return 0;
}

/* Parse "f:a,f:a,..." into tones. Returns count, or -1 on bad input. */
static int parse_tones(const char *text, tone_t *tones, size_t max_tones) {
    if (text == NULL) {
        return 0;
    }
    size_t count = 0;
    const char *p = text;
    while (*p != '\0') {
        char *end = NULL;
        const double freq = strtod(p, &end);
        if (end == p || *end != ':') {
            return -1;
        }
        p = end + 1;
        const double amp = strtod(p, &end);
        if (end == p) {
            return -1;
        }
        p = end;
        if (count >= max_tones) {
            return -1;
        }
        tones[count].freq = freq;
        tones[count].amp = amp;
        tones[count].phase = 0.0;
        count++;
        if (*p == ',') {
            p++;
        } else if (*p != '\0') {
            return -1;
        }
    }
    return (int)count;
}

static int cmd_synthesize(int argc, char **argv) {
    const char *type = opt_value(argc, argv, "type");
    const char *output = opt_value(argc, argv, "output");
    if (type == NULL || output == NULL) {
        fprintf(stderr,
                "fourier-lab synthesize: --type and --output required\n");
        return 2;
    }
    double fs = 0.0;
    if (opt_double(argc, argv, "fs", &fs) < 0 || fs <= 0.0) {
        fprintf(stderr, "fourier-lab synthesize: valid --fs required\n");
        return 2;
    }
    double duration = 0.0;
    if (opt_double(argc, argv, "duration", &duration) < 0 ||
        duration <= 0.0) {
        fprintf(stderr,
                "fourier-lab synthesize: valid --duration required\n");
        return 2;
    }
    const size_t n = (size_t)(fs * duration + 0.5);
    if (n == 0) {
        fprintf(stderr, "fourier-lab synthesize: zero samples requested\n");
        return 2;
    }
    double *x = (double *)malloc(n * sizeof(double));
    if (x == NULL) {
        fprintf(stderr, "fourier-lab synthesize: out of memory\n");
        return 1;
    }
    int rc = 0;
    if (strcmp(type, "tones") == 0) {
        tone_t tones[64];
        const int nt = parse_tones(opt_value(argc, argv, "tones"), tones, 64);
        if (nt < 0) {
            fprintf(stderr,
                    "fourier-lab synthesize: bad --tones, expected "
                    "\"freq:amp[,freq:amp...]\"\n");
            free(x);
            return 2;
        }
        if (nt == 0) {
            fprintf(stderr, "fourier-lab synthesize: --tones required for "
                            "type 'tones'\n");
            free(x);
            return 2;
        }
        rc = signal_tones(x, n, fs, tones, (size_t)nt);
        printf("synthesize: %d tone(s), fs=%.17g Hz, %.6g s, %zu samples\n",
               nt, fs, duration, n);
    } else if (strcmp(type, "square") == 0 || strcmp(type, "sawtooth") == 0) {
        double freq = 0.0;
        double amp = 1.0;
        long harmonics = 25;
        if (opt_double(argc, argv, "freq", &freq) < 0 || freq <= 0.0) {
            fprintf(stderr, "fourier-lab synthesize: valid --freq required\n");
            free(x);
            return 2;
        }
        if (opt_double(argc, argv, "amp", &amp) < 0) {
            fprintf(stderr, "fourier-lab synthesize: invalid --amp value\n");
            free(x);
            return 2;
        }
        if (opt_long(argc, argv, "harmonics", &harmonics) < 0 ||
            harmonics < 1) {
            fprintf(stderr,
                    "fourier-lab synthesize: invalid --harmonics value\n");
            free(x);
            return 2;
        }
        if (strcmp(type, "square") == 0) {
            rc = signal_square(x, n, fs, freq, amp, (int)harmonics);
            printf("synthesize: square wave %.6g Hz, amp %.6g, %ld odd "
                   "harmonics, fs=%.17g Hz, %zu samples\n",
                   freq, amp, harmonics, fs, n);
        } else {
            rc = signal_sawtooth(x, n, fs, freq, amp, (int)harmonics);
            printf("synthesize: sawtooth wave %.6g Hz, amp %.6g, %ld "
                   "harmonics, fs=%.17g Hz, %zu samples\n",
                   freq, amp, harmonics, fs, n);
        }
    } else if (strcmp(type, "ecg") == 0) {
        double bpm = 72.0;
        double noise_std = 0.0;
        long seed = 42;
        if (opt_double(argc, argv, "bpm", &bpm) < 0 || bpm <= 0.0) {
            fprintf(stderr, "fourier-lab synthesize: invalid --bpm value\n");
            free(x);
            return 2;
        }
        if (opt_double(argc, argv, "noise-std", &noise_std) < 0 ||
            noise_std < 0.0) {
            fprintf(stderr,
                    "fourier-lab synthesize: invalid --noise-std value\n");
            free(x);
            return 2;
        }
        if (opt_long(argc, argv, "seed", &seed) < 0) {
            fprintf(stderr, "fourier-lab synthesize: invalid --seed value\n");
            free(x);
            return 2;
        }
        rc = signal_ecg(x, n, fs, bpm, noise_std, (uint64_t)seed);
        printf("synthesize: ECG %.6g BPM, noise std %.6g, seed %ld, "
               "fs=%.17g Hz, %zu samples\n",
               bpm, noise_std, seed, fs, n);
    } else if (strcmp(type, "noise") == 0) {
        double std = 1.0;
        long seed = 42;
        if (opt_double(argc, argv, "noise-std", &std) < 0 || std < 0.0) {
            fprintf(stderr,
                    "fourier-lab synthesize: invalid --noise-std value\n");
            free(x);
            return 2;
        }
        if (opt_long(argc, argv, "seed", &seed) < 0) {
            fprintf(stderr, "fourier-lab synthesize: invalid --seed value\n");
            free(x);
            return 2;
        }
        rc = signal_noise(x, n, std, (uint64_t)seed);
        printf("synthesize: white Gaussian noise, std %.6g, seed %ld, "
               "fs=%.17g Hz, %zu samples\n",
               std, seed, fs, n);
    } else if (strcmp(type, "chirp") == 0) {
        double f0 = 0.0;
        double f1 = 0.0;
        double amp = 1.0;
        if (opt_double(argc, argv, "f0", &f0) < 0 || f0 < 0.0) {
            fprintf(stderr, "fourier-lab synthesize: valid --f0 required\n");
            free(x);
            return 2;
        }
        if (opt_double(argc, argv, "f1", &f1) < 0 || f1 < 0.0) {
            fprintf(stderr, "fourier-lab synthesize: valid --f1 required\n");
            free(x);
            return 2;
        }
        if (opt_double(argc, argv, "amp", &amp) < 0) {
            fprintf(stderr, "fourier-lab synthesize: invalid --amp value\n");
            free(x);
            return 2;
        }
        rc = signal_chirp(x, n, fs, f0, f1, amp);
        printf("synthesize: chirp %.6g -> %.6g Hz, amp %.6g, fs=%.17g Hz, "
               "%zu samples\n",
               f0, f1, amp, fs, n);
    } else {
        fprintf(stderr, "fourier-lab synthesize: unknown type '%s'\n", type);
        free(x);
        return 2;
    }
    if (rc != 0) {
        fprintf(stderr, "fourier-lab synthesize: generation failed\n");
        free(x);
        return 1;
    }
    if (io_write_signal_csv(output, x, n, fs) != 0) {
        fprintf(stderr, "fourier-lab synthesize: cannot write '%s'\n", output);
        free(x);
        return 1;
    }
    printf("  written to %s\n", output);
    free(x);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        print_usage(stderr);
        return 2;
    }
    const char *cmd = argv[1];
    if (strcmp(cmd, "--help") == 0 || strcmp(cmd, "-h") == 0 ||
        strcmp(cmd, "help") == 0) {
        print_usage(stdout);
        return 0;
    }
    if (strcmp(cmd, "--version") == 0 || strcmp(cmd, "-V") == 0) {
        printf("fourier-lab %s\n", FOURIER_LAB_VERSION);
        return 0;
    }
    /* Subcommand sees only its own arguments (argv[2..]). */
    argc -= 2;
    argv += 2;
    if (strcmp(cmd, "fft") == 0) {
        return cmd_fft(argc, argv);
    }
    if (strcmp(cmd, "analyze") == 0) {
        return cmd_analyze(argc, argv);
    }
    if (strcmp(cmd, "synthesize") == 0) {
        return cmd_synthesize(argc, argv);
    }
    if (strcmp(cmd, "reconstruct") == 0) {
        return cmd_reconstruct(argc, argv);
    }
    if (strcmp(cmd, "qrs") == 0) {
        return cmd_qrs(argc, argv);
    }
    fprintf(stderr, "fourier-lab: unknown command '%s'\n", cmd);
    print_usage(stderr);
    return 2;
}
