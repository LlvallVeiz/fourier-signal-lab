/* test_main.c -- Test runner: runs every suite and exits non-zero
 * when any check fails.
 */

#include <stdio.h>

#include "test_framework.h"

/* Single definition of the shared check counters. */
int fl_checks = 0;
int fl_failures = 0;

void test_fft(void);
void test_window(void);
void test_signal(void);
void test_spectrum(void);
void test_io(void);

int main(void) {
    int rc = 0;

    test_fft();
    rc |= fl_report("test_fft");
    fl_checks = 0;
    fl_failures = 0;

    test_window();
    rc |= fl_report("test_window");
    fl_checks = 0;
    fl_failures = 0;

    test_signal();
    rc |= fl_report("test_signal");
    fl_checks = 0;
    fl_failures = 0;

    test_spectrum();
    rc |= fl_report("test_spectrum");
    fl_checks = 0;
    fl_failures = 0;

    test_io();
    rc |= fl_report("test_io");

    if (rc == 0) {
        printf("ALL TESTS PASSED\n");
    } else {
        printf("SOME TESTS FAILED\n");
    }
    return rc;
}
