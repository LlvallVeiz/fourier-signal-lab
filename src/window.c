/* window.c -- Data windows (rectangular, Hann, Hamming, Blackman).
 *
 * Periodic (DFT) definitions, w[n] for n = 0..N-1:
 *   rectangular: 1
 *   Hann:       0.5 - 0.5*cos(2*pi*n / N)
 *   Hamming:    0.54 - 0.46*cos(2*pi*n / N)
 *   Blackman:   0.42 - 0.5*cos(2*pi*n / N) + 0.08*cos(4*pi*n / N)
 * See Harris (1978) for the trade-offs and tabulated figures of merit.
 */

#include "window.h"

#include <ctype.h>
#include <math.h>
#include <string.h>

#include "complex.h"

int window_parse(const char *name, window_t *out) {
    if (name == NULL || out == NULL) {
        return -1;
    }
    /* Lower-case copy for a case-insensitive match. */
    char buf[32];
    size_t i = 0;
    while (name[i] != '\0' && i < sizeof(buf) - 1) {
        buf[i] = (char)tolower((unsigned char)name[i]);
        i++;
    }
    buf[i] = '\0';
    if (strcmp(buf, "rect") == 0 || strcmp(buf, "rectangular") == 0 ||
        strcmp(buf, "none") == 0) {
        *out = WINDOW_RECTANGULAR;
    } else if (strcmp(buf, "hann") == 0 || strcmp(buf, "hanning") == 0) {
        *out = WINDOW_HANN;
    } else if (strcmp(buf, "hamming") == 0) {
        *out = WINDOW_HAMMING;
    } else if (strcmp(buf, "blackman") == 0) {
        *out = WINDOW_BLACKMAN;
    } else {
        return -1;
    }
    return 0;
}

const char *window_name(window_t w) {
    switch (w) {
        case WINDOW_RECTANGULAR: return "rectangular";
        case WINDOW_HANN:        return "hann";
        case WINDOW_HAMMING:     return "hamming";
        case WINDOW_BLACKMAN:    return "blackman";
        default:                 return "unknown";
    }
}

double window_coef(window_t w, size_t n, size_t N) {
    if (N == 0 || n >= N) {
        return 0.0;
    }
    const double phase = 2.0 * FL_PI * (double)n / (double)N;
    switch (w) {
        case WINDOW_RECTANGULAR:
            return 1.0;
        case WINDOW_HANN:
            return 0.5 - 0.5 * cos(phase);
        case WINDOW_HAMMING:
            return 0.54 - 0.46 * cos(phase);
        case WINDOW_BLACKMAN:
            return 0.42 - 0.5 * cos(phase) + 0.08 * cos(2.0 * phase);
        default:
            return 0.0;
    }
}

void window_apply(window_t w, double *x, size_t n) {
    if (x == NULL) {
        return;
    }
    for (size_t i = 0; i < n; i++) {
        x[i] *= window_coef(w, i, n);
    }
}

double window_coherent_gain(window_t w, size_t n) {
    if (n == 0) {
        return 0.0;
    }
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += window_coef(w, i, n);
    }
    return sum / (double)n;
}

double window_enbw(window_t w, size_t n) {
    if (n == 0) {
        return 0.0;
    }
    double sum = 0.0;
    double sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        const double c = window_coef(w, i, n);
        sum += c;
        sum_sq += c * c;
    }
    if (sum == 0.0) {
        return 0.0;
    }
    return (double)n * sum_sq / (sum * sum);
}

double window_amplitude_correction(window_t w, size_t n) {
    const double gain = window_coherent_gain(w, n);
    if (gain == 0.0) {
        return 0.0;
    }
    return 1.0 / gain;
}

double window_highest_sidelobe_db(window_t w) {
    switch (w) {
        case WINDOW_RECTANGULAR: return -13.0;
        case WINDOW_HANN:        return -31.0;
        case WINDOW_HAMMING:     return -43.0;
        case WINDOW_BLACKMAN:    return -58.0;
        default:                 return 0.0;
    }
}
