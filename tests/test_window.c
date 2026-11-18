/* test_window.c -- Tests for window coefficients and figures of merit.
 *
 * Reference values (Harris, 1978, Table I; also obtained analytically
 * for the periodic definitions used here):
 *   window      coherent gain   ENBW (bins)   highest sidelobe
 *   rectangular      1.00            1.00           -13 dB
 *   Hann             0.50            1.50           -31 dB
 *   Hamming          0.54            1.36           -43 dB
 *   Blackman         0.42            1.73           -58 dB
 * The coherent gain and ENBW are checked against the analytic values
 * (exact for the periodic window: e.g. Hann sum = N/2, so ENBW = 1.5
 * exactly; Hamming ENBW = 0.3974/0.2916 = 1.3628; Blackman ENBW =
 * 0.3046/0.1764 = 1.7268).
 */

#include "test_framework.h"

#include <math.h>

#include "window.h"

void test_window(void) {
    const size_t n = 1024;

    /* --- Rectangular: all ones, gain 1, ENBW 1, no correction. --- */
    {
        for (size_t i = 0; i < n; i++) {
            CHECK_NEAR(window_coef(WINDOW_RECTANGULAR, i, n), 1.0, 1e-15);
        }
        CHECK_NEAR(window_coherent_gain(WINDOW_RECTANGULAR, n), 1.0, 1e-12);
        CHECK_NEAR(window_enbw(WINDOW_RECTANGULAR, n), 1.0, 1e-12);
        CHECK_NEAR(window_amplitude_correction(WINDOW_RECTANGULAR, n), 1.0,
                   1e-12);
        CHECK_NEAR(window_highest_sidelobe_db(WINDOW_RECTANGULAR), -13.0,
                   1e-9);
    }

    /* --- Hann: sum N/2, endpoints 0, centre 1, ENBW 1.50. --- */
    {
        CHECK_NEAR(window_coef(WINDOW_HANN, 0, n), 0.0, 1e-15);
        CHECK_NEAR(window_coef(WINDOW_HANN, n / 2, n), 1.0, 1e-12);
        CHECK_NEAR(window_coef(WINDOW_HANN, n / 4, n), 0.5, 1e-12);
        double sum = 0.0;
        for (size_t i = 0; i < n; i++) {
            sum += window_coef(WINDOW_HANN, i, n);
        }
        CHECK_NEAR(sum, (double)n / 2.0, 1e-9);
        CHECK_NEAR(window_coherent_gain(WINDOW_HANN, n), 0.5, 1e-12);
        CHECK_NEAR(window_enbw(WINDOW_HANN, n), 1.5, 1e-9);
        CHECK_NEAR(window_amplitude_correction(WINDOW_HANN, n), 2.0, 1e-12);
        CHECK_NEAR(window_highest_sidelobe_db(WINDOW_HANN), -31.0, 1e-9);
    }

    /* --- Hamming: endpoints 0.08 (the pedestal), ENBW ~1.3628. --- */
    {
        CHECK_NEAR(window_coef(WINDOW_HAMMING, 0, n), 0.08, 1e-12);
        CHECK_NEAR(window_coef(WINDOW_HAMMING, n / 2, n), 1.0, 1e-12);
        CHECK_NEAR(window_coherent_gain(WINDOW_HAMMING, n), 0.54, 1e-12);
        CHECK_NEAR(window_enbw(WINDOW_HAMMING, n), 1.3628, 0.001);
        CHECK_NEAR(window_amplitude_correction(WINDOW_HAMMING, n),
                   1.0 / 0.54, 1e-9);
        CHECK_NEAR(window_highest_sidelobe_db(WINDOW_HAMMING), -43.0, 1e-9);
    }

    /* --- Blackman: endpoints 0, centre 1, ENBW ~1.7268. --- */
    {
        CHECK_NEAR(window_coef(WINDOW_BLACKMAN, 0, n), 0.0, 1e-12);
        CHECK_NEAR(window_coef(WINDOW_BLACKMAN, n / 2, n), 1.0, 1e-12);
        CHECK_NEAR(window_coherent_gain(WINDOW_BLACKMAN, n), 0.42, 1e-12);
        CHECK_NEAR(window_enbw(WINDOW_BLACKMAN, n), 1.7268, 0.001);
        CHECK_NEAR(window_amplitude_correction(WINDOW_BLACKMAN, n),
                   1.0 / 0.42, 1e-9);
        CHECK_NEAR(window_highest_sidelobe_db(WINDOW_BLACKMAN), -58.0, 1e-9);
    }

    /* --- window_apply scales a constant signal by the window sum. --- */
    {
        double x[16];
        for (int i = 0; i < 16; i++) {
            x[i] = 1.0;
        }
        window_apply(WINDOW_HANN, x, 16);
        double sum = 0.0;
        for (int i = 0; i < 16; i++) {
            sum += x[i];
        }
        CHECK_NEAR(sum, 8.0, 1e-12);
        CHECK_NEAR(x[0], 0.0, 1e-15);
        CHECK_NEAR(x[8], 1.0, 1e-12);
    }

    /* --- Parsing: names, aliases, case, and rejection. --- */
    {
        window_t w;
        CHECK(window_parse("hann", &w) == 0 && w == WINDOW_HANN);
        CHECK(window_parse("Hanning", &w) == 0 && w == WINDOW_HANN);
        CHECK(window_parse("HAMMING", &w) == 0 && w == WINDOW_HAMMING);
        CHECK(window_parse("blackman", &w) == 0 && w == WINDOW_BLACKMAN);
        CHECK(window_parse("rect", &w) == 0 && w == WINDOW_RECTANGULAR);
        CHECK(window_parse("rectangular", &w) == 0 &&
              w == WINDOW_RECTANGULAR);
        CHECK(window_parse("none", &w) == 0 && w == WINDOW_RECTANGULAR);
        CHECK(window_parse("kaiser", &w) == -1);
        CHECK(window_parse("", &w) == -1);
        CHECK(window_parse(NULL, &w) == -1);
        CHECK(window_parse("hann", NULL) == -1);
    }

    /* --- Names round-trip through the parser. --- */
    {
        for (int wi = 0; wi < WINDOW_COUNT; wi++) {
            const window_t w = (window_t)wi;
            window_t back;
            CHECK(window_parse(window_name(w), &back) == 0 && back == w);
        }
    }

    /* --- Out-of-range coefficient requests are safe. --- */
    CHECK_NEAR(window_coef(WINDOW_HANN, 5, 0), 0.0, 1e-15);
    CHECK_NEAR(window_coef(WINDOW_HANN, 16, 16), 0.0, 1e-15);
    CHECK_NEAR(window_coherent_gain(WINDOW_HANN, 0), 0.0, 1e-15);
    CHECK_NEAR(window_enbw(WINDOW_HANN, 0), 0.0, 1e-15);
}
