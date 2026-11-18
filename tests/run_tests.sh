#!/usr/bin/env bash
# run_tests.sh -- CLI smoke tests for the fourier-lab binary.
# Assumes the binary is already built (make) and run_tests passed.
set -euo pipefail
cd "$(dirname "$0")/.."

./fourier-lab --version
./fourier-lab --help > /dev/null

# synthesize -> fft -> analyze -> reconstruct round trip
./fourier-lab synthesize --type tones --tones "50:1.0,120:0.5" \
    --fs 1024 --duration 1.0 --output /tmp/fl_smoke_tones.csv
./fourier-lab fft --input /tmp/fl_smoke_tones.csv \
    --output /tmp/fl_smoke_spec.csv --window hann
./fourier-lab analyze --input /tmp/fl_smoke_tones.csv --window hann --top 5
./fourier-lab reconstruct --input /tmp/fl_smoke_tones.csv
./fourier-lab synthesize --type square --freq 10 --amp 1.0 --harmonics 25 \
    --fs 1024 --duration 1.0 --output /tmp/fl_smoke_square.csv
./fourier-lab analyze --input /tmp/fl_smoke_square.csv --window hann --top 5
./fourier-lab synthesize --type sawtooth --freq 10 --amp 1.0 --harmonics 50 \
    --fs 1024 --duration 1.0 --output /tmp/fl_smoke_saw.csv
./fourier-lab synthesize --type ecg --bpm 72 --fs 360 --duration 10 \
    --output /tmp/fl_smoke_ecg.csv
./fourier-lab qrs --input /tmp/fl_smoke_ecg.csv
./fourier-lab synthesize --type noise --noise-std 1.0 --seed 7 \
    --fs 1000 --duration 1.0 --output /tmp/fl_smoke_noise.csv
./fourier-lab synthesize --type chirp --f0 10 --f1 100 --amp 1.0 \
    --fs 1000 --duration 1.0 --output /tmp/fl_smoke_chirp.csv

# Expected failures: unknown command, unknown window, missing file,
# non-power-of-two N, bad tones syntax.
if ./fourier-lab frobnicate 2>/dev/null; then
    echo "smoke: unknown command should have failed" >&2; exit 1
fi
if ./fourier-lab fft --input /tmp/fl_smoke_tones.csv \
        --output /tmp/fl_smoke_x.csv --window kaiser 2>/dev/null; then
    echo "smoke: unknown window should have failed" >&2; exit 1
fi
if ./fourier-lab analyze --input /tmp/fl_no_such_file.csv 2>/dev/null; then
    echo "smoke: missing file should have failed" >&2; exit 1
fi
if ./fourier-lab fft --input /tmp/fl_smoke_tones.csv \
        --output /tmp/fl_smoke_x.csv --n 100 2>/dev/null; then
    echo "smoke: non-power-of-two N should have failed" >&2; exit 1
fi
if ./fourier-lab synthesize --type tones --tones "50x1.0" --fs 1024 \
        --duration 1.0 --output /tmp/fl_smoke_x.csv 2>/dev/null; then
    echo "smoke: bad tones syntax should have failed" >&2; exit 1
fi

rm -f /tmp/fl_smoke_*.csv /tmp/fl_smoke_x.csv
echo "CLI smoke tests passed."
