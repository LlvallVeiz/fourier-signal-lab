#!/usr/bin/env python3
"""validate_numpy.py -- Cross-check the C FFT against numpy.fft.

Generates a two-tone signal with the C CLI, computes its spectrum with
the C CLI (rectangular and Hann windows), and compares magnitudes and
phases against numpy.fft.rfft on the same samples. Exits non-zero when
the worst relative difference exceeds 1e-9 (rectangular) or 1e-6 (Hann,
whose window coefficients differ in the last ulp between the two
implementations).

Requires numpy and the built ./fourier-lab binary (run from the repo
root, or set FOURIER_LAB_BIN).
"""

import math
import os
import subprocess
import sys
import tempfile

import numpy as np

BIN = os.environ.get("FOURIER_LAB_BIN", "./fourier-lab")
FS = 1024.0
N = 1024


def run(*args):
    result = subprocess.run([BIN, *args], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"FAILED: {' '.join(args)}\n{result.stderr}")
        sys.exit(1)
    return result.stdout


def read_signal(path):
    xs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            xs.append(float(line))
    return np.array(xs)


def read_spectrum(path):
    freq, mag, phase = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("freq"):
                continue
            a, b, c = line.split(",")
            freq.append(float(a))
            mag.append(float(b))
            phase.append(float(c))
    return np.array(freq), np.array(mag), np.array(phase)


def check_window(tmp, window, tol):
    sig = os.path.join(tmp, f"tones_{window}.csv")
    spec = os.path.join(tmp, f"spec_{window}.csv")
    run("synthesize", "--type", "tones", "--tones", "50:1.0,120:0.5,311:0.25",
        "--fs", str(int(FS)), "--duration", "1.0", "--output", sig)
    run("fft", "--input", sig, "--output", spec, "--window", window)

    x = read_signal(sig)
    assert len(x) == N, f"expected {N} samples, got {len(x)}"
    freq, mag, phase = read_spectrum(spec)

    # Reference: numpy rfft of the same windowed samples.
    if window == "rectangular":
        w = np.ones(N)
        corr = 1.0
    elif window == "hann":
        w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(N) / N)
        corr = 2.0
    else:  # pragma: no cover - only rect/hann are checked here
        raise ValueError(window)
    ref = np.fft.rfft(x * w)
    ref_mag = np.abs(ref) * corr * 2.0 / N
    ref_mag[0] /= 2.0
    ref_mag[-1] /= 2.0
    ref_phase = np.angle(ref)

    mag_err = np.max(np.abs(mag - ref_mag))
    # Phase is only meaningful where the magnitude is non-negligible.
    strong = ref_mag > 1e-6
    phase_err = np.max(np.abs(
        (phase[strong] - ref_phase[strong] + np.pi) % (2 * np.pi) - np.pi))
    print(f"window={window}: max |mag diff| = {mag_err:.3e}, "
          f"max |phase diff| (strong bins) = {phase_err:.3e}")
    assert mag_err < tol, f"magnitude mismatch for {window}: {mag_err}"
    assert phase_err < 1e-9, f"phase mismatch for {window}: {phase_err}"

    # The three tones must read back at their true amplitudes.
    for f0, amp in ((50.0, 1.0), (120.0, 0.5), (311.0, 0.25)):
        k = int(round(f0))
        got = mag[k]
        assert math.isclose(got, amp, rel_tol=1e-9, abs_tol=1e-9), \
            f"tone {f0} Hz: got {got}, expected {amp}"
    print(f"window={window}: tone amplitudes 50 Hz=1.0, 120 Hz=0.5, "
          f"311 Hz=0.25 -- OK")


def main():
    if not os.path.isfile(BIN):
        print(f"binary not found: {BIN} (run make first)")
        sys.exit(1)
    with tempfile.TemporaryDirectory() as tmp:
        check_window(tmp, "rectangular", 1e-9)
        check_window(tmp, "hann", 1e-6)
    print("numpy validation PASSED")


if __name__ == "__main__":
    main()
