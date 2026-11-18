/* test_framework.h -- Minimal assert-based test framework (no external
 * dependencies). Each check increments a counter; failures print the
 * expression, file and line. Suites register their checks through
 * CHECK / CHECK_NEAR and report a per-suite total.
 */

#ifndef FOURIER_LAB_TEST_FRAMEWORK_H
#define FOURIER_LAB_TEST_FRAMEWORK_H

#include <math.h>
#include <stdio.h>

/* Counters live in exactly one translation unit (test_main.c); every
 * suite shares them through these extern declarations. (A `static`
 * counter in this header would silently give each test file its own
 * copy and every suite would report 0/0.) */
extern int fl_checks;
extern int fl_failures;

static inline void fl_check(int cond, const char *expr, const char *file,
                            int line) {
    fl_checks++;
    if (!cond) {
        fl_failures++;
        fprintf(stderr, "FAIL %s:%d: %s\n", file, line, expr);
    }
}

static inline void fl_check_near(double actual, double expected, double tol,
                                 const char *expr, const char *file,
                                 int line) {
    fl_checks++;
    if (!(fabs(actual - expected) <= tol)) {
        fl_failures++;
        fprintf(stderr, "FAIL %s:%d: %s -- got %.17g, expected %.17g "
                        "(+/- %.3g)\n",
                file, line, expr, actual, expected, tol);
    }
}

static inline int fl_report(const char *suite) {
    printf("%s: %d/%d checks passed", suite, fl_checks - fl_failures,
           fl_checks);
    if (fl_failures == 0) {
        printf(" -- OK\n");
        return 0;
    }
    printf(" -- %d FAILURES\n", fl_failures);
    return 1;
}

#define CHECK(cond) fl_check((cond), #cond, __FILE__, __LINE__)
#define CHECK_NEAR(actual, expected, tol) \
    fl_check_near((actual), (expected), (tol), \
                  #actual " ~= " #expected, __FILE__, __LINE__)

#endif /* FOURIER_LAB_TEST_FRAMEWORK_H */
