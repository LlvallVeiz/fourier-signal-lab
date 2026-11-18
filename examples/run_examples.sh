#!/usr/bin/env bash
# run_examples.sh -- Regenerate the example outputs under examples/outputs/.
# Every number printed below comes from a real run of the CLI; the
# headline figures are collected in examples/RESULTS.md.
set -euo pipefail
cd "$(dirname "$0")/.."
BIN=./fourier-lab
OUT=examples/outputs
mkdir -p "$OUT"

echo "=== 1. Two tones: 50 Hz (amp 1.0) + 120 Hz (amp 0.5) ==="
$BIN synthesize --type tones --tones "50:1.0,120:0.5" \
    --fs 1024 --duration 1.0 --output "$OUT/two_tones.csv"
$BIN analyze --input "$OUT/two_tones.csv" --window hann --top 5

echo
echo "=== 2. Synthetic ECG, 72 BPM ==="
$BIN synthesize --type ecg --bpm 72 --fs 360 --duration 10 \
    --output "$OUT/ecg_72bpm.csv"
$BIN qrs --input "$OUT/ecg_72bpm.csv"
$BIN analyze --input "$OUT/ecg_72bpm.csv" --window hann --top 5

echo
echo "=== 3. Square wave, 10 Hz, 25 odd harmonics ==="
$BIN synthesize --type square --freq 10 --amp 1.0 --harmonics 25 \
    --fs 1024 --duration 1.0 --output "$OUT/square_10hz.csv"
$BIN analyze --input "$OUT/square_10hz.csv" --window hann --top 5

echo
echo "=== 4. Spectral leakage: off-bin tone, 55.3 Hz ==="
$BIN synthesize --type tones --tones "55.3:1.0" \
    --fs 1024 --duration 1.0 --output "$OUT/offbin_55hz.csv"
echo "--- rectangular window ---"
$BIN analyze --input "$OUT/offbin_55hz.csv" --window rectangular --top 3
echo "--- Hann window ---"
$BIN analyze --input "$OUT/offbin_55hz.csv" --window hann --top 3

echo
echo "=== 5. FFT -> IFFT reconstruction (two tones) ==="
$BIN reconstruct --input "$OUT/two_tones.csv"

echo
echo "=== 6. Aliasing: 700 Hz tone sampled at 1024 Hz ==="
$BIN synthesize --type tones --tones "700:1.0" \
    --fs 1024 --duration 1.0 --output "$OUT/alias_700hz.csv"
$BIN analyze --input "$OUT/alias_700hz.csv" --window rectangular --top 3

echo
echo "=== 7. Timing: FFT of 2^20 samples ==="
$BIN synthesize --type noise --noise-std 1.0 --seed 7 \
    --fs 1048576 --duration 1.0 --output "$OUT/noise_1m.csv"
time $BIN fft --input "$OUT/noise_1m.csv" \
    --output "$OUT/noise_1m_spec.csv" --window hann
# The 2^20-sample files are ~46 MB of text; their numbers are recorded
# in RESULTS.md, so drop them again to keep the tree lean.
rm -f "$OUT/noise_1m.csv" "$OUT/noise_1m_spec.csv"

echo
echo "=== 8. Chirp sweeping 10 -> 100 Hz ==="
$BIN synthesize --type chirp --f0 10 --f1 100 --amp 1.0 \
    --fs 1000 --duration 1.0 --output "$OUT/chirp.csv"
$BIN analyze --input "$OUT/chirp.csv" --window hann --top 5

echo
echo "Examples done. Outputs in $OUT/"
