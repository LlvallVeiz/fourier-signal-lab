# Example Results

All figures below are the actual output of `bash examples/run_examples.sh`
(run 2026-11-17) plus `make bench` and `python3 tests/validate_numpy.py`.
Nothing here is hand-computed or estimated, except where a value is
explicitly labelled as a textbook prediction.

## 1. Two tones: 50 Hz (amp 1.0) + 120 Hz (amp 0.5)

`synthesize --type tones --tones "50:1.0,120:0.5" --fs 1024 --duration 1.0`,
then `analyze --window hann`:

| Peak | Frequency | Magnitude |
|------|-----------|-----------|
| 1 | 50.000 Hz | 1.000 |
| 2 | 120.000 Hz | 0.500 |

Both tones are recovered at exactly their true amplitudes (the Hann
amplitude correction of 2.0 does its job). Parseval check:
sum(x²) = (1/N)·sum(|X|²) = 240, relative difference 1.37e-14.

Note: the reported "leakage" of 0.2000 for this signal is **not** window
leakage — it is the second tone's energy, which lies outside ±2 bins of
the dominant peak by definition. The leakage metric is only meaningful
for a single tone (see §4).

## 2. Synthetic ECG, 72 BPM

`synthesize --type ecg --bpm 72 --fs 360 --duration 10`, then `qrs`:

- R peaks detected: **12** (ground truth: 12 beats in 10 s at 72 BPM)
- R times (s): 0.4167, 1.2500, 2.0833, 2.9167, 3.7500, 4.5833, 5.4167,
  6.2500, 7.0833, 7.9167, 8.7500, 9.5833 — exactly the beat positions
  the generator placed (0.4167 + k·0.8333 s)
- Heart rate: **72.00 BPM** (ground truth 72)

Spectral analysis of the same signal (Hann, N=4096, df=0.0879 Hz) finds
peaks at harmonics of the 1.2 Hz beat rate — 2.402 Hz (0.0730),
3.600 Hz (0.0839), 6.002 Hz (0.0633), 7.199 Hz (0.0728) — plus a small
near-DC peak at 0.044 Hz (0.0830) from the signal's non-zero mean. The
lesson: for heart rate, the time-domain QRS detector is the right tool;
the spectrum spreads the pulse train's energy over many harmonics.

## 3. Square wave, 10 Hz, 25 odd harmonics

`synthesize --type square --freq 10 --amp 1.0 --harmonics 25`, then
`analyze --window hann`. The odd harmonics follow the textbook 4/(πk)
law (amplitude 1.0):

| Harmonic | Frequency | Measured | Theory 4/(πk) |
|----------|-----------|----------|---------------|
| 1st | 10 Hz | 1.27324 | 1.27324 |
| 3rd | 30 Hz | 0.424413 | 0.424413 |
| 5th | 50 Hz | 0.254648 | 0.254648 |
| 7th | 70 Hz | 0.181891 | 0.181891 |
| 9th | 90 Hz | 0.141471 | 0.141471 |

Even harmonics (20, 40 Hz) read < 0.01, as theory requires. Parseval
relative difference: 1.21e-14.

## 4. Spectral leakage: off-bin tone, 55.3 Hz

A tone at 55.3 Hz does not complete an integer number of cycles in the
1024-sample record, so it leaks. Same signal, two windows:

| Window | Detected frequency | Magnitude | Leakage (energy outside ±2 bins) |
|--------|--------------------|-----------|----------------------------------|
| rectangular | 55.134 Hz (error 0.166 Hz) | 0.858 | 0.0531 |
| Hann | 55.316 Hz (error 0.016 Hz) | 0.943 | 0.0002 |

The Hann window cuts the leaked energy by a factor of ~230 and improves
the parabolically interpolated frequency estimate tenfold, at the cost
of a wider main lobe (ENBW 1.50 bins vs 1.00) — the classic trade-off
tabulated by Harris (1978).

## 5. FFT → IFFT reconstruction

`reconstruct` on the two-tone signal (N=1024):

- Max absolute error: **2.73e-14**
- RMS error: **9.06e-15**
- Relative max error (vs max |x| = 1.4933): **1.83e-14**

The round trip is exact to double-precision rounding.

## 6. Aliasing: 700 Hz tone sampled at 1024 Hz

Nyquist is 512 Hz, so 700 Hz cannot be represented. The analyser
reports a single peak at **324.000 Hz** (1024 − 700) with magnitude
1.000 — the alias. There is no way to tell from the samples alone
whether the tone was 324 Hz or 700 Hz: that is the Nyquist–Shannon
sampling theorem in action.

## 7. Timing

`make bench` (this machine, two runs):

- FFT, N = 2²⁰ (1,048,576 samples): **≈ 0.07 s**
- DFT, N = 1024: **0.031 s**; FFT, N = 1024: **0.000029 s** — the FFT is
  **~1100× faster** at N = 1024 (1075×–1182× across runs)
- Extrapolating the DFT's O(N²) cost to N = 2²⁰ gives ~9 hours; the FFT
  does it in 0.07 s

The `fft` subcommand on the 2²⁰-sample noise file took 12.8 s of wall
time in the example run — that is dominated by reading/writing the
~1M-line / ~524k-line text CSV files, not by the transform (0.07 s).
Text I/O, not the FFT, is the bottleneck of the CLI at this size.

## 8. Chirp sweeping 10 → 100 Hz

`synthesize --type chirp --f0 10 --f1 100 --fs 1000 --duration 1.0`:
the spectrum is spread over the whole 10–100 Hz band, so the analyser
finds a single broad peak at 56.244 Hz (magnitude 0.206) and reports
leakage 0.8595 — correctly, since a chirp is not a tone. A single FFT
of a chirp answers "which frequencies are present", not "when".

## Cross-check against numpy

`python3 tests/validate_numpy.py` (numpy 2.5.4) compares the C FFT
against `numpy.fft.rfft` on the same three-tone signal
(50 Hz / 1.0, 120 Hz / 0.5, 311 Hz / 0.25):

| Window | Max \|magnitude difference\| | Max \|phase difference\| |
|--------|---------------------------|------------------------|
| rectangular | 8.2e-15 | 0.0 |
| Hann | 8.1e-15 | 4.4e-16 |

All three tone amplitudes read back exactly (1.0 / 0.5 / 0.25).
