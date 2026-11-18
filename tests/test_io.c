/* test_io.c -- Tests for the signal/spectrum CSV readers and writers.
 * Round trips go through /tmp so the repo tree stays clean.
 */

#include "test_framework.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "io.h"

void test_io(void) {
    const char *sig_path = "/tmp/fourier_lab_test_signal.csv";
    const char *spec_path = "/tmp/fourier_lab_test_spectrum.csv";

    /* --- Signal round trip preserves values and fs exactly. --- */
    {
        double x[5] = {0.0, 1.5, -2.25, 3.141592653589793, 1e-12};
        CHECK(io_write_signal_csv(sig_path, x, 5, 1000.0) == 0);
        double *y = NULL;
        size_t n = 0;
        double fs = 0.0;
        CHECK(io_read_signal_csv(sig_path, &y, &n, &fs) == 0);
        CHECK(n == 5);
        CHECK_NEAR(fs, 1000.0, 1e-9);
        for (int i = 0; i < 5; i++) {
            CHECK_NEAR(y[i], x[i], 1e-15);
        }
        free(y);
    }

    /* --- A file without an fs header reports fs = -1. --- */
    {
        FILE *f = fopen(sig_path, "w");
        CHECK(f != NULL);
        fprintf(f, "1.0\n2.0\n3.0\n");
        fclose(f);
        double *y = NULL;
        size_t n = 0;
        double fs = 0.0;
        CHECK(io_read_signal_csv(sig_path, &y, &n, &fs) == 0);
        CHECK(n == 3);
        CHECK_NEAR(fs, -1.0, 1e-9);
        CHECK_NEAR(y[1], 2.0, 1e-15);
        free(y);
    }

    /* --- Comments and blank lines are skipped. --- */
    {
        FILE *f = fopen(sig_path, "w");
        CHECK(f != NULL);
        fprintf(f, "# fs=44100\n# a comment\n\n1.0\n\n2.0\n");
        fclose(f);
        double *y = NULL;
        size_t n = 0;
        double fs = 0.0;
        CHECK(io_read_signal_csv(sig_path, &y, &n, &fs) == 0);
        CHECK(n == 2);
        CHECK_NEAR(fs, 44100.0, 1e-9);
        free(y);
    }

    /* --- Malformed data and missing files fail cleanly. --- */
    {
        FILE *f = fopen(sig_path, "w");
        CHECK(f != NULL);
        fprintf(f, "1.0\nhello\n");
        fclose(f);
        double *y = NULL;
        size_t n = 0;
        double fs = 0.0;
        CHECK(io_read_signal_csv(sig_path, &y, &n, &fs) == -1);
        CHECK(io_read_signal_csv("/tmp/fourier_lab_no_such_file.csv", &y,
                                 &n, &fs) == -1);
        CHECK(io_read_signal_csv(NULL, &y, &n, &fs) == -1);
        CHECK(io_write_signal_csv("/tmp/no_such_dir/x.csv", y, 0, 1.0) == -1);
    }

    /* --- Spectrum CSV writes a header and one row per bin. --- */
    {
        double freq[3] = {0.0, 1.0, 2.0};
        double mag[3] = {0.0, 2.0, 0.5};
        double phase[3] = {0.0, -1.5707963267948966, 0.0};
        CHECK(io_write_spectrum_csv(spec_path, freq, mag, phase, 3, 1024.0,
                                    1024, "hann") == 0);
        FILE *f = fopen(spec_path, "r");
        CHECK(f != NULL);
        char line[256];
        CHECK(fgets(line, sizeof(line), f) != NULL);
        CHECK(fgets(line, sizeof(line), f) != NULL); /* column header */
        double fr = 0.0, mg = 0.0, ph = 0.0;
        CHECK(fscanf(f, "%lf,%lf,%lf", &fr, &mg, &ph) == 3);
        CHECK_NEAR(fr, 0.0, 1e-12);
        CHECK(fscanf(f, "%lf,%lf,%lf", &fr, &mg, &ph) == 3);
        CHECK_NEAR(fr, 1.0, 1e-12);
        CHECK_NEAR(mg, 2.0, 1e-12);
        CHECK_NEAR(ph, -1.5707963267948966, 1e-12);
        fclose(f);
        CHECK(io_write_spectrum_csv(spec_path, NULL, mag, phase, 3, 1024.0,
                                    1024, "hann") == -1);
    }

    remove(sig_path);
    remove(spec_path);
}
