/* io.c -- Plain-text signal and spectrum files. */

#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int io_read_signal_csv(const char *path, double **out_x, size_t *out_n,
                       double *out_fs) {
    if (path == NULL || out_x == NULL || out_n == NULL || out_fs == NULL) {
        return -1;
    }
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }
    size_t cap = 1024;
    size_t n = 0;
    double *x = (double *)malloc(cap * sizeof(double));
    if (x == NULL) {
        fclose(f);
        return -1;
    }
    *out_fs = -1.0;
    char line[256];
    int ok = 0;
    while (fgets(line, sizeof(line), f) != NULL) {
        if (line[0] == '#') {
            double fs = 0.0;
            if (sscanf(line, "# fs=%lf", &fs) == 1) {
                *out_fs = fs;
            }
            continue;
        }
        /* Skip blank lines. */
        char *p = line;
        while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') {
            p++;
        }
        if (*p == '\0') {
            continue;
        }
        char *end = NULL;
        const double v = strtod(line, &end);
        if (end == line) {
            goto fail; /* malformed line */
        }
        if (n == cap) {
            cap *= 2;
            double *nx = (double *)realloc(x, cap * sizeof(double));
            if (nx == NULL) {
                goto fail;
            }
            x = nx;
        }
        x[n++] = v;
    }
    ok = 1;
fail:
    fclose(f);
    if (!ok) {
        free(x);
        return -1;
    }
    *out_x = x;
    *out_n = n;
    return 0;
}

int io_write_signal_csv(const char *path, const double *x, size_t n,
                        double fs) {
    if (path == NULL || x == NULL) {
        return -1;
    }
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        return -1;
    }
    fprintf(f, "# fs=%.17g\n", fs);
    for (size_t i = 0; i < n; i++) {
        fprintf(f, "%.17g\n", x[i]);
    }
    if (fclose(f) != 0) {
        return -1;
    }
    return 0;
}

int io_write_spectrum_csv(const char *path, const double *freq,
                          const double *mag, const double *phase,
                          size_t nbins, double fs, size_t n,
                          const char *window) {
    if (path == NULL || freq == NULL || mag == NULL || phase == NULL ||
        window == NULL) {
        return -1;
    }
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        return -1;
    }
    fprintf(f, "# fs=%.17g, n=%zu, window=%s\n", fs, n, window);
    fprintf(f, "freq,magnitude,phase\n");
    for (size_t k = 0; k < nbins; k++) {
        fprintf(f, "%.17g,%.17g,%.17g\n", freq[k], mag[k], phase[k]);
    }
    if (fclose(f) != 0) {
        return -1;
    }
    return 0;
}
