/* spectrum.c -- Spectrum reading, peak detection, Parseval, leakage,
 * and a simplified Pan-Tompkins QRS detector.
 */

#include "spectrum.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void spectrum_magnitude(const cx_t *X, size_t n, window_t w, double *mag) {
    if (X == NULL || mag == NULL || n == 0) {
        return;
    }
    const double corr = window_amplitude_correction(w, n);
    const size_t nbins = n / 2 + 1;
    for (size_t k = 0; k < nbins; k++) {
        double scale = 2.0 / (double)n;
        if (k == 0 || k == n / 2) {
            scale = 1.0 / (double)n; /* DC / Nyquist: no mirror image */
        }
        mag[k] = scale * corr * cx_abs(X[k]);
    }
}

void spectrum_phase(const cx_t *X, size_t n, double *phase) {
    if (X == NULL || phase == NULL || n == 0) {
        return;
    }
    const size_t nbins = n / 2 + 1;
    for (size_t k = 0; k < nbins; k++) {
        phase[k] = cx_arg(X[k]);
    }
}

/* Insertion-sort helper: sort candidate indices by magnitude, largest
 * first. Candidates are bin indices stored in idx[]. */
static void sort_by_magnitude_desc(const double *mag, size_t *idx,
                                   size_t count) {
    for (size_t i = 1; i < count; i++) {
        const size_t key = idx[i];
        size_t j = i;
        while (j > 0 && mag[idx[j - 1]] < mag[key]) {
            idx[j] = idx[j - 1];
            j--;
        }
        idx[j] = key;
    }
}

/* Insertion-sort helper: sort accepted peaks by bin ascending. */
static void sort_peaks_by_bin(peak_t *peaks, size_t count) {
    for (size_t i = 1; i < count; i++) {
        const peak_t key = peaks[i];
        size_t j = i;
        while (j > 0 && peaks[j - 1].bin > key.bin) {
            peaks[j] = peaks[j - 1];
            j--;
        }
        peaks[j] = key;
    }
}

int spectrum_find_peaks(const double *mag, size_t nbins, double df,
                        double min_height, size_t min_distance,
                        peak_t *peaks, size_t max_peaks) {
    if (mag == NULL || peaks == NULL || nbins < 3 || max_peaks == 0) {
        return -1;
    }
    /* Collect local maxima above the height threshold (skip DC). */
    size_t *cand = (size_t *)malloc(nbins * sizeof(size_t));
    if (cand == NULL) {
        return -1;
    }
    size_t ncand = 0;
    for (size_t k = 1; k + 1 < nbins; k++) {
        if (mag[k] >= min_height && mag[k] >= mag[k - 1] &&
            mag[k] > mag[k + 1]) {
            cand[ncand++] = k;
        }
    }
    sort_by_magnitude_desc(mag, cand, ncand);

    size_t nfound = 0;
    for (size_t i = 0; i < ncand && nfound < max_peaks; i++) {
        const size_t k = cand[i];
        int too_close = 0;
        for (size_t j = 0; j < nfound; j++) {
            const size_t d = peaks[j].bin > k ? peaks[j].bin - k : k - peaks[j].bin;
            if (d < min_distance) {
                too_close = 1;
                break;
            }
        }
        if (too_close) {
            continue;
        }
        peaks[nfound].bin = k;
        peaks[nfound].freq = spectrum_parabolic_freq(mag, k, df);
        peaks[nfound].magnitude = mag[k];
        peaks[nfound].phase = 0.0; /* filled in by the caller if needed */
        nfound++;
    }
    free(cand);
    sort_peaks_by_bin(peaks, nfound);
    return (int)nfound;
}

double spectrum_parabolic_delta(const double *mag, size_t k) {
    double a, b, c;
    if (mag[k - 1] > 0.0 && mag[k] > 0.0 && mag[k + 1] > 0.0) {
        /* Log-magnitude parabola (see spectrum.h). */
        a = log(mag[k - 1]);
        b = log(mag[k]);
        c = log(mag[k + 1]);
    } else {
        a = mag[k - 1];
        b = mag[k];
        c = mag[k + 1];
    }
    const double denom = a - 2.0 * b + c;
    if (fabs(denom) < 1e-300) {
        return 0.0;
    }
    double delta = 0.5 * (a - c) / denom;
    if (delta > 1.0) {
        delta = 1.0;
    } else if (delta < -1.0) {
        delta = -1.0;
    }
    return delta;
}

double spectrum_parabolic_freq(const double *mag, size_t k, double df) {
    return ((double)k + spectrum_parabolic_delta(mag, k)) * df;
}

double spectrum_parseval_time(const double *x, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

double spectrum_parseval_freq(const cx_t *X, size_t n) {
    double sum = 0.0;
    for (size_t k = 0; k < n; k++) {
        sum += X[k].re * X[k].re + X[k].im * X[k].im;
    }
    return sum / (double)n;
}

double spectrum_leakage(const cx_t *X, size_t n, size_t spread) {
    if (X == NULL || n == 0) {
        return 0.0;
    }
    /* Dominant bin over the single-sided spectrum. */
    size_t peak = 0;
    double peak_e = 0.0;
    double total = 0.0;
    for (size_t k = 0; k <= n / 2; k++) {
        const double e = X[k].re * X[k].re + X[k].im * X[k].im;
        total += e;
        if (e > peak_e) {
            peak_e = e;
            peak = k;
        }
    }
    if (total == 0.0) {
        return 0.0;
    }
    double main_e = 0.0;
    const size_t lo = peak > spread ? peak - spread : 0;
    const size_t hi = peak + spread < n / 2 ? peak + spread : n / 2;
    for (size_t k = lo; k <= hi; k++) {
        main_e += X[k].re * X[k].re + X[k].im * X[k].im;
    }
    return 1.0 - main_e / total;
}

int spectrum_qrs_detect(const double *x, size_t n, double fs,
                        size_t *r_peaks, size_t max_peaks) {
    if (x == NULL || r_peaks == NULL || n < 16 || fs <= 0.0 ||
        max_peaks == 0) {
        return -1;
    }
    double *mwi = (double *)malloc(n * sizeof(double));
    double *sq = (double *)malloc(n * sizeof(double));
    if (mwi == NULL || sq == NULL) {
        free(mwi);
        free(sq);
        return -1;
    }
    /* 1. Five-point derivative (Pan-Tompkins): emphasises the steep
     * QRS slope, suppresses the slow P/T waves. */
    for (size_t i = 0; i < n; i++) {
        if (i < 2 || i + 2 >= n) {
            sq[i] = 0.0;
        } else {
            const double d =
                (2.0 * x[i + 2] + x[i + 1] - x[i - 1] - 2.0 * x[i - 2]) / 8.0;
            sq[i] = d * d; /* 2. Squaring: all positive, QRS dominates. */
        }
    }
    /* 3. Moving-window integration over 150 ms. */
    const size_t win = (size_t)(0.15 * fs + 0.5);
    const size_t half = win / 2;
    double max_mwi = 0.0;
    for (size_t i = 0; i < n; i++) {
        const size_t lo = i > half ? i - half : 0;
        const size_t hi = i + half < n ? i + half : n - 1;
        double sum = 0.0;
        for (size_t j = lo; j <= hi; j++) {
            sum += sq[j];
        }
        mwi[i] = sum / (double)(hi - lo + 1);
        if (mwi[i] > max_mwi) {
            max_mwi = mwi[i];
        }
    }
    free(sq);
    if (max_mwi == 0.0) {
        free(mwi);
        return 0;
    }
    /* 4. Threshold + refractory period (200 ms). The moving-window
     * integral of a QRS is a flat-topped hump (the window is wider
     * than the complex), so the first local maximum sits on the
     * hump's leading edge, ~half a window early. Refine each hit by
     * taking the largest raw sample within +-100 ms: the R spike
     * dominates everything else in that neighbourhood, so the refined
     * index is the true R position. */
    const double threshold = 0.25 * max_mwi;
    const size_t refractory = (size_t)(0.20 * fs + 0.5);
    const size_t refine = (size_t)(0.10 * fs + 0.5);
    size_t count = 0;
    size_t last = 0;
    int have_last = 0;
    for (size_t i = 1; i + 1 < n; i++) {
        if (mwi[i] < threshold || mwi[i] < mwi[i - 1] ||
            mwi[i] < mwi[i + 1]) {
            continue; /* not a local maximum above threshold */
        }
        if (have_last && i - last < refractory) {
            continue;
        }
        if (count >= max_peaks) {
            break;
        }
        const size_t lo = i > refine ? i - refine : 0;
        const size_t hi = i + refine < n ? i + refine : n - 1;
        size_t best = lo;
        for (size_t j = lo + 1; j <= hi; j++) {
            if (x[j] > x[best]) {
                best = j;
            }
        }
        r_peaks[count++] = best;
        last = i;
        have_last = 1;
    }
    free(mwi);
    return (int)count;
}

double spectrum_bpm(const size_t *r_peaks, size_t count, double fs) {
    if (r_peaks == NULL || count < 2 || fs <= 0.0) {
        return -1.0;
    }
    const double mean_rr =
        (double)(r_peaks[count - 1] - r_peaks[0]) / (double)(count - 1) / fs;
    if (mean_rr <= 0.0) {
        return -1.0;
    }
    return 60.0 / mean_rr;
}
