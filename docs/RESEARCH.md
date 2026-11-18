# Research Notes — fourier-signal-lab

How the Fourier transform, windowing, sampling, and QRS detection
work, why this lab implements them the way it does, and where every
claim comes from. All references below were read (abstract/full text)
before writing; no citation here is invented.

## 1. The discrete Fourier transform and the FFT

For a length-N sampled signal x[n], the discrete Fourier transform is

    X[k] = sum_{n=0}^{N-1} x[n] * exp(-2*pi*i*k*n / N),  k = 0..N-1

computed directly, that is N² complex multiply-adds. Cooley and Tukey
(1965) showed that when N is composite, the DFT splits into smaller
DFTs of the even- and odd-indexed samples, and recursing the split
gives O(N log N) work — for N = 2²⁰, about fifty thousand times fewer
operations than the direct sum (N² = 2⁴⁰ vs N·log2 N = 20·2²⁰). Their
paper reports a 2048-point complex DFT in 0.02 minutes on an IBM 7094.
The same idea had been used by Gauss around 1805 to interpolate
asteroid orbits, but it was Cooley and Tukey's paper that made it a
practical tool.

This lab implements the iterative radix-2 decimation-in-time form
(`src/fft.c`): a bit-reversal permutation followed by log2(N) butterfly
stages, in place, O(1) extra storage. A naive O(N²) DFT
(`dft_forward`) is kept as a reference implementation, and the test
suite checks the FFT against it on random real and complex inputs
(max difference < 1e-9 at N=256).

Measured on this machine (`make bench`): FFT at N = 2²⁰ takes ~0.07 s;
at N = 1024 the FFT is ~1075× faster than the naive DFT (0.000029 s vs
0.031 s). Extrapolating the DFT's quadratic cost to N = 2²⁰ gives
~9 hours — the same order of magnitude as the gap Cooley and Tukey
closed in 1965.

The inverse transform is implemented by conjugation
(IFFT(x) = conj(FFT(conj(x))) / N), and the round trip FFT → IFFT
recovers the input to 2.7e-14 max absolute error (see
`examples/RESULTS.md`).

## 2. Windows, leakage, and figures of merit

A DFT implicitly treats the finite record as one period of an infinite
periodic signal. A tone that does not complete an integer number of
cycles in the record is discontinuous at the boundary, and its energy
smears across bins — spectral leakage. Multiplying the record by a
window that tapers to zero at the ends suppresses the discontinuity at
the cost of a wider main lobe (worse frequency resolution) and a
coherent-gain loss that must be corrected when reading amplitudes.

Harris (1978) is the standard catalogue: dozens of windows with their
sidelobe levels, equivalent noise bandwidth (ENBW), scalloping loss,
and worst-case processing loss tabulated side by side. This lab
implements the four classic windows and computes the figures of merit
numerically from the coefficients, so they are exact for the N used:

| Window | Definition (periodic) | Coherent gain | ENBW (bins) | Highest sidelobe |
|--------|----------------------|---------------|-------------|------------------|
| rectangular | 1 | 1.00 | 1.00 | −13 dB |
| Hann | 0.5 − 0.5·cos(2πn/N) | 0.50 | 1.50 | −31 dB |
| Hamming | 0.54 − 0.46·cos(2πn/N) | 0.54 | 1.36 | −43 dB |
| Blackman | 0.42 − 0.5·cos(2πn/N) + 0.08·cos(4πn/N) | 0.42 | 1.73 | −58 dB |

Coherent gain and ENBW are verified in the test suite against the
analytic values (e.g. Hann ENBW = 1.50 exactly for the periodic
definition); the sidelobe figures are Harris's Table I values
(rounded). The amplitude correction is 1/coherent gain (2.0 for Hann),
and the test suite confirms that an on-bin tone reads its true
amplitude under every window.

Measured on a 55.3 Hz off-bin tone (N=1024): rectangular leakage
(energy outside ±2 bins of the peak) is 0.0531, Hann is 0.0002 — a
factor of ~230 less. The interpolated frequency estimate improves from
55.134 Hz (error 0.166 Hz) to 55.316 Hz (error 0.016 Hz).

## 3. Peak interpolation

A bin is only accurate to ±df/2. Fitting a parabola through the three
bins around a magnitude peak refines the estimate. Smith ("Spectral
Audio Signal Processing", CCRMA, Stanford) shows that fitting the
parabola to the log (dB) magnitude is about twice as accurate as
fitting linear magnitude, and exact for a Gaussian window (whose
transform is a parabola on a dB scale). This lab fits log-magnitude
(`spectrum_parabolic_delta`), with a fallback to linear magnitude when
a bin is zero. Measured: a 50.3 Hz tone is estimated at 50.316 Hz with
a Hann window (error 0.016 Hz = 0.016 bins) and at 50.134 Hz with a
rectangular window (error 0.166 Hz) — the rectangular window's
Dirichlet-kernel main lobe is a poor parabola, which is why the
rectangular tolerance in the tests is loose.

## 4. Fourier series: square and sawtooth waves

A square wave of amplitude A is the sum of its odd harmonics only:

    x(t) = (4A/π) · Σ_{odd k} sin(2πkft) / k

and a sawtooth is the sum of all harmonics with alternating signs:

    x(t) = (2A/π) · Σ_{k≥1} (−1)^(k+1) sin(2πkft) / k

The lab synthesises both from these series (`signal_square`,
`signal_sawtooth`), so the harmonic amplitudes follow the 1/k law
exactly — verified in the examples: a 10 Hz square wave reads 1.27324
(= 4/π), 0.424413 (= 4/3π), 0.254648 (= 4/5π) at 10/30/50 Hz, with even
harmonics below 0.01.

Truncating the series at a jump leaves the Gibbs phenomenon: the
partial sums overshoot the jump by ~8.95% of its height (the
Wilbraham–Gibbs constant), no matter how many harmonics are kept —
first computed by Wilbraham (1848), rediscovered by Gibbs (1899).
Measured here: a 25-harmonic square wave of amplitude 1.0 peaks at
1.1791, i.e. 0.179 above 1.0 for a jump of height 2 — 8.95% of the
jump, matching the constant to three decimals.

## 5. Sampling and aliasing

Nyquist (1928) showed that a channel of bandwidth B carries at most
2B independent pulses per second; Shannon (1949, "Communication in
the Presence of Noise", Theorem 1) proved that a signal containing no
frequencies above W is completely determined by samples spaced
1/(2W) apart. A signal sampled below that rate folds: a 700 Hz tone
sampled at 1024 Hz is indistinguishable from a 324 Hz tone
(1024 − 700). The lab demonstrates this directly (`examples/RESULTS.md`
§6): the analyser reports a single peak at 324.000 Hz with magnitude
1.000. Nothing in the samples can tell the two apart — which is why
real acquisition chains put an anti-aliasing filter before the ADC.

## 6. Synthetic ECG and the Pan–Tompkins detector

Pan and Tompkins (1985) detect QRS complexes in real time with a
pipeline of band-pass filtering, differentiation, squaring, moving
window integration, and adaptive thresholds; on the MIT/BIH arrhythmia
database it correctly detects 99.3% of QRS complexes. The lab's
synthetic ECG is a teaching model in the same spirit: each heartbeat
is a sum of Gaussian bumps — P (0.15 mV), Q (−0.10 mV), R (1.00 mV),
S (−0.25 mV), T (0.30 mV) — with a narrow QRS (~25 ms) on a slow
baseline. It reproduces the morphology a QRS detector keys on and
nothing more; it is not a physiological simulation.

The detector (`spectrum_qrs_detect`) is the textbook pipeline reduced
to its essentials: a 5-point derivative, squaring, 150 ms moving
window integration, a threshold at 25% of the maximum, a 200 ms
refractory period, and a refinement step that takes the largest raw
sample within ±100 ms of each hit (the R spike). On the synthetic ECG
(72 BPM, 10 s) it finds all 12 beats at exactly the generated
positions (0.4167 + k·0.8333 s) and reports 72.00 BPM — including with
0.05 mV of added noise. The full Pan–Tompkins adaptive thresholds and
T-wave discrimination are deliberately not implemented; the README
lists them under limitations.

The ECG spectrum itself is a useful negative lesson: the pulse train's
energy spreads over harmonics of the 1.2 Hz beat rate (peaks at 2.4,
3.6, 6.0, 7.2 Hz), so the spectrum is the wrong tool for heart rate —
the time-domain detector is the right one.

## 7. Validation strategy

Three independent checks guard against a wrong transform:

1. **Hand values**: DFT of [1,0,0,0] is [1,1,1,1]; DFT of [0,1,0,0]
   is exp(−2πik/4); a sine exactly on a bin has |X[k]| = A·N/2;
   Parseval's identity sum(x²) = (1/N)·sum(|X|²) holds to 1e-12.
2. **Cross-implementation**: the FFT agrees with the naive DFT to
   < 1e-9 (N=256, complex input).
3. **External reference**: `tests/validate_numpy.py` compares the C
   FFT against `numpy.fft.rfft` on the same samples — max magnitude
   difference 8.2e-15 (rectangular) and 8.1e-15 (Hann), phase
   difference ≤ 4.4e-16.

## References

1. Cooley, J. W., and Tukey, J. W. "An Algorithm for the Machine
   Calculation of Complex Fourier Series." *Mathematics of
   Computation* 19(90), pp. 297–301, 1965.
   DOI: 10.1090/S0025-5718-1965-0178586-1.
   https://doi.org/10.1090/S0025-5718-1965-0178586-1
2. Harris, F. J. "On the Use of Windows for Harmonic Analysis with
   the Discrete Fourier Transform." *Proceedings of the IEEE* 66(1),
   pp. 51–83, 1978. DOI: 10.1109/PROC.1978.10837.
   Full text (PDF): http://web.mit.edu/~xiphmont/Public/windows.pdf
3. Pan, J., and Tompkins, W. J. "A Real-Time QRS Detection
   Algorithm." *IEEE Transactions on Biomedical Engineering* BME-32(3),
   pp. 230–236, 1985. DOI: 10.1109/TBME.1985.325532. PMID: 3997178.
   Full text (PDF):
   https://people.ece.cornell.edu/land/courses/ece5030/labs/s2014/Pan_Tompkins.pdf
4. Shannon, C. E. "Communication in the Presence of Noise."
   *Proceedings of the Institute of Radio Engineers* 37(1), pp. 10–21,
   1949. DOI: 10.1109/jrproc.1949.232969.
   Full text (PDF):
   https://www.math.unl.edu/~bdeng1/Teaching/math943/Topics/Genetics/Shannon1949.pdf
5. Nyquist, H. "Certain Topics in Telegraph Transmission Theory."
   *Bell System Technical Journal* 7(2), pp. 280–308, 1928.
6. Smith, J. O., III. *Spectral Audio Signal Processing*. W3K
   Publishing, 2011. ISBN 978-0-9745607-3-1. Online edition, CCRMA,
   Stanford University — chapter "Quadratic Interpolation of Spectral
   Peaks":
   https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
7. Wilbraham, H. "On a Certain Periodic Function." *Cambridge and
   Dublin Mathematical Journal* 3, pp. 198–201, 1848; Gibbs, J. W.
   "Fourier's Series." *Nature* 59, p. 606, 1899. (Gibbs phenomenon:
   ~8.95% overshoot at a jump discontinuity.)
   Overview: https://en.wikipedia.org/wiki/Gibbs_phenomenon
8. NumPy Developers. `numpy.fft` — used only as an external
   cross-check in `tests/validate_numpy.py`; the lab itself has no
   dependencies.
