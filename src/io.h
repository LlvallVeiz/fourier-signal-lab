/* io.h -- Plain-text signal and spectrum files.
 *
 * Signal CSV: one sample per line, plain decimal. An optional comment
 * header "# fs=<rate>" records the sampling rate; comment lines
 * starting with '#' are ignored on read.
 *
 * Spectrum CSV: header "# fs=<rate>, n=<N>, window=<name>", then one
 * "freq,magnitude,phase" row per bin.
 *
 * All functions return 0 on success, -1 on error (file not found,
 * malformed data, or allocation failure).
 */

#ifndef FOURIER_LAB_IO_H
#define FOURIER_LAB_IO_H

#include <stddef.h>

/* Read a signal CSV. Allocates *out_x (caller frees with free()).
 * *out_n is the sample count; *out_fs is the rate from the "# fs="
 * header, or -1.0 when the file has no header. */
int io_read_signal_csv(const char *path, double **out_x, size_t *out_n,
                       double *out_fs);

/* Write a signal CSV with a "# fs=" header. */
int io_write_signal_csv(const char *path, const double *x, size_t n,
                        double fs);

/* Write a spectrum CSV: freq (Hz), magnitude, phase (rad). */
int io_write_spectrum_csv(const char *path, const double *freq,
                          const double *mag, const double *phase,
                          size_t nbins, double fs, size_t n,
                          const char *window);

#endif /* FOURIER_LAB_IO_H */
