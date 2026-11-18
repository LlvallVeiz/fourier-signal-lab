# fourier-signal-lab

Laboratorium analisis sinyal dalam C11 murni, tanpa dependensi
eksternal: FFT radix-2 Cooley–Tukey dari nol, DFT naif sebagai
referensi, jendela spektral, sintesis sinyal (nada, gelombang kotak /
gigi gergaji dari deret Fourier, ECG sintetik, derau, chirp), analisis
spektral (deteksi puncak + interpolasi parabolik), dan detektor QRS
ala Pan–Tompkins yang disederhanakan.

## Latar belakang

Transformasi Fourier diskret (DFT) menguraikan sinyal menjadi
kandungan frekuensinya:

    X[k] = Σ x[n] · exp(−2πi·k·n / N)

Dihitung langsung, biayanya N² operasi. Cooley dan Tukey (1965)
menunjukkan bahwa DFT dapat dipecah menjadi DFT-DFT yang lebih kecil
dari sampel genap dan ganjil, dan dengan mengulang pemecahan itu
biayanya menjadi O(N log N) — untuk N = 2²⁰, sekitar lima puluh ribu
kali lebih sedikit operasi daripada penjumlahan langsung. Lab ini
mengimplementasikan bentuk iteratif radix-2 (bit-reversal + butterfly
in-place) dan menyimpan DFT naif O(N²) sebagai referensi yang
divalidasi silang oleh tes.

Dua efek yang wajib dipahami siapa pun yang memakai FFT ikut
didemonstrasikan di sini:

- **Kebocoran spektral (leakage).** DFT menganggap rekaman berhingga
  sebagai satu periode sinyal periodik. Nada yang tidak menyelesaikan
  jumlah siklus bulat di dalam rekaman menjadi diskontinu di batas,
  dan energinya menyebar ke bin tetangga. Mengalikan rekaman dengan
  jendela (Hann, Hamming, Blackman) menekan diskontinuitas itu dengan
  harga main lobe yang lebih lebar — trade-off yang dikatalogkan
  Harris (1978).
- **Aliasing.** Shannon (1949) membuktikan sinyal tanpa frekuensi di
  atas W ditentukan penuh oleh sampel berjarak 1/(2W). Nada 700 Hz yang
  dicuplik pada 1024 Hz tidak dapat dibedakan dari nada 324 Hz
  (1024 − 700) — dan memang begitulah hasil yang dilaporkan lab ini.

Detail teori, rumus, dan verifikasi hitungan tangan ada di
[`docs/RESEARCH.md`](docs/RESEARCH.md). Semua angka contoh di README
ini berasal dari run beneran (`examples/RESULTS.md`).

## Struktur proyek

```
fourier-signal-lab/
├── src/
│   ├── complex.h            # aritmetika kompleks minimal (static inline)
│   ├── fft.h / fft.c        # DFT naif + FFT radix-2 iteratif + IFFT
│   ├── window.h / window.c  # rectangular, Hann, Hamming, Blackman +
│   │                        # coherent gain / ENBW / koreksi amplitudo
│   ├── signal.h / signal.c  # sintesis: nada, kotak, gigi gergaji,
│   │                        # ECG sintetik, derau, chirp (+ PRNG seed)
│   ├── spectrum.h / spectrum.c  # magnitudo/fase, deteksi puncak,
│   │                        # interpolasi parabolik (log-magnitudo),
│   │                        # Parseval, leakage, detektor QRS
│   ├── io.h / io.c          # file CSV sinyal & spektrum
│   └── main.c               # CLI: fft | analyze | synthesize |
│                            #   reconstruct | qrs
├── tests/
│   ├── test_framework.h     # framework assert ringan (tanpa dependensi)
│   ├── test_fft.c           # DFT hitung tangan, FFT vs DFT, Parseval
│   ├── test_window.c        # koefisien & angka Harris (ENBW, gain)
│   ├── test_signal.c        # nilai hitung tangan: sinus, kotak (Gibbs),
│   │                        # gigi gergaji, ECG, derau, chirp
│   ├── test_spectrum.c      # puncak dua-nada, interpolasi, leakage,
│   │                        # aliasing, harmonik kotak, deteksi QRS
│   ├── test_io.c            # round-trip CSV
│   ├── test_main.c          # runner
│   ├── run_tests.sh         # smoke test CLI
│   └── validate_numpy.py    # FFT C vs numpy.fft.rfft
├── examples/
│   ├── run_examples.sh      # menjalankan semua contoh
│   ├── benchmark_fft.c      # benchmark FFT vs DFT (`make bench`)
│   ├── RESULTS.md           # angka hasil run beneran
│   └── outputs/             # hasil generate (diabaikan git)
├── docs/
│   └── RESEARCH.md          # catatan riset + References
├── Makefile                 # make / make test / make examples / make bench
├── CMakeLists.txt           # alternatif build + ctest
├── LICENSE                  # MIT
└── README.md
```

## Cara jalan

```sh
make            # build binary ./fourier-lab (gcc, -Wall -Wextra, 0 warning)
make test       # tes unit + validasi numpy + smoke test CLI
make examples   # jalankan semua contoh ke examples/outputs/
make bench      # benchmark FFT vs DFT
```

Atau dengan CMake:

```sh
cmake -S . -B build && cmake --build build
cd build && ctest
```

Kebutuhan: compiler C11 (gcc/clang) dan `make`. Validasi numpy butuh
Python 3 + numpy (opsional — tes C tetap jalan tanpanya; `ctest` akan
melewati tes numpy bila Python tidak ada).

## Penggunaan CLI

Sinyal adalah file CSV teks: satu sampel per baris, dengan header
opsional `# fs=<rate>`.

```sh
# Sintesis sinyal dua-nada, lalu analisis spektrumnya
./fourier-lab synthesize --type tones --tones "50:1.0,120:0.5" \
    --fs 1024 --duration 1.0 --output two_tones.csv
./fourier-lab analyze --input two_tones.csv --window hann --top 5

# Tulis spektrum ke CSV (freq, magnitude, phase)
./fourier-lab fft --input two_tones.csv --output spec.csv --window hann

# Rekonstruksi FFT -> IFFT dan laporkan errornya
./fourier-lab reconstruct --input two_tones.csv

# ECG sintetik 72 BPM + deteksi QRS
./fourier-lab synthesize --type ecg --bpm 72 --fs 360 --duration 10 \
    --output ecg.csv
./fourier-lab qrs --input ecg.csv
```

Jenis sintesis: `tones` (`--tones "f:a,f:a"`), `square` / `sawtooth`
(`--freq --amp --harmonics`), `ecg` (`--bpm --noise-std --seed`),
`noise` (`--noise-std --seed`), `chirp` (`--f0 --f1 --amp`).
Jendela: `rectangular`, `hann`, `hamming`, `blackman` (default `hann`).

## Contoh output (dari run beneran)

### Dua nada: 50 Hz (amp 1,0) + 120 Hz (amp 0,5)

```
analyze: examples/outputs/two_tones.csv
  N=1024  fs=1024 Hz  df=1 Hz  Nyquist=512 Hz  window=hann
  window: coherent gain 0.5000, ENBW 1.5000 bins, amplitude correction 2.0000, highest sidelobe -31 dB
  peaks (freq, magnitude):
        50.000 Hz  1
       120.000 Hz  0.5
  Parseval: sum(x^2)=240  (1/N)sum(|X|^2)=240  rel. diff=1.37e-14
```

Kedua nada terbaca tepat pada amplitudo aslinya — koreksi amplitudo
Hann (×2,0) bekerja. Cek Parseval: jumlah kuadrat di domain waktu dan
frekuensi sepakat sampai 1,37e-14.

### ECG sintetik 72 BPM + detektor QRS

```
qrs: examples/outputs/ecg_72bpm.csv (fs=360 Hz, 3600 samples)
  R peaks detected: 12
  R times (s): 0.4167 1.2500 2.0833 2.9167 3.7500 4.5833 5.4167 6.2500 7.0833 7.9167 8.7500 9.5833
  heart rate: 72.00 BPM
```

12 detak terdeteksi, semuanya tepat di posisi yang digenerate
(0,4167 + k·0,8333 dtk), heart rate 72,00 BPM — sama dengan ground
truth. Spektrum sinyal yang sama justru menyebar ke harmonik laju
detak (2,4 / 3,6 / 6,0 / 7,2 Hz): untuk heart rate, detektor QRS
domain-waktu adalah alat yang benar, bukan spektrum.

### Gelombang kotak 10 Hz: harmonik mengikuti hukum 1/k

| Harmonik | Frekuensi | Terukur | Teori 4/(πk) |
|----------|-----------|---------|--------------|
| 1 | 10 Hz | 1,27324 | 1,27324 |
| 3 | 30 Hz | 0,424413 | 0,424413 |
| 5 | 50 Hz | 0,254648 | 0,254648 |
| 7 | 70 Hz | 0,181891 | 0,181891 |
| 9 | 90 Hz | 0,141471 | 0,141471 |

Harmonik genap terbaca < 0,01, sesuai teori. Overshoot Gibbs pada
gelombang kotak 25-harmonik terukur 1,1791 — 8,95% dari lompatan
setinggi 2, persis konstanta Wilbraham–Gibbs.

### Kebocoran spektral: nada 55,3 Hz (di luar bin)

| Jendela | Frekuensi terdeteksi | Magnitudo | Leakage (energi di luar ±2 bin) |
|---------|----------------------|-----------|---------------------------------|
| rectangular | 55,134 Hz (galat 0,166 Hz) | 0,858 | 0,0531 |
| Hann | 55,316 Hz (galat 0,016 Hz) | 0,943 | 0,0002 |

Jendela Hann memangkas energi bocor ~230× dan memperbaiki estimasi
frekuensi sepuluh kali lipat, dengan harga main lobe yang lebih lebar
(ENBW 1,50 bin vs 1,00).

### Rekonstruksi FFT → IFFT

```
reconstruct: examples/outputs/two_tones.csv (N=1024, fs=1024 Hz)
  FFT -> IFFT round trip over 1024 samples
  max abs error: 2.73115e-14
  RMS error:     9.05811e-15
  max |x|:       1.4933  (relative max error 1.82893e-14)
```

### Aliasing: nada 700 Hz dicuplik pada 1024 Hz

Analyser melaporkan satu puncak pada **324,000 Hz** (1024 − 700)
dengan magnitudo 1,000. Dari sampelnya saja, tidak ada cara membedakan
nada 324 Hz dan 700 Hz.

### Kecepatan

`make bench` (mesin ini, dua run): FFT N = 2²⁰ (1.048.576 sampel)
**≈ 0,07 dtk**; pada N = 1024, FFT 0,000029 dtk vs DFT 0,031 dtk —
**~1100× lebih cepat** (1075–1182× antar-run). Ekstrapolasi biaya
O(N²) DFT ke N = 2²⁰ menghasilkan ~9 jam. (Waktu wall 12,8 dtk pada
perintah `fft` untuk file 2²⁰ sampel didominasi I/O CSV teks, bukan
transformasinya.)

### Validasi vs numpy

`tests/validate_numpy.py` (numpy 2.5.4) membandingkan FFT C dengan
`numpy.fft.rfft` pada sinyal tiga-nada yang sama:

| Jendela | Selisih magnitudo maks | Selisih fase maks |
|---------|------------------------|-------------------|
| rectangular | 8,2e-15 | 0,0 |
| Hann | 8,1e-15 | 4,4e-16 |

Ketiga amplitudo nada terbaca tepat (1,0 / 0,5 / 0,25).

## Verifikasi

```
test_fft:      130/130 checks passed -- OK
test_window:   1070/1070 checks passed -- OK
test_signal:   1228/1228 checks passed -- OK
test_spectrum: 68/68 checks passed -- OK
test_io:       34/34 checks passed -- OK
ALL TESTS PASSED
numpy validation PASSED
CLI smoke tests passed.
```

Tes mencakup nilai hitung tangan (DFT impuls `[1,0,0,0]` → semua 1;
sinus tepat di bin → |X[k]| = A·N/2; Parseval; rekonstruksi ~1e-14;
jumlah koefisien jendela Hann = N/2; ENBW vs tabel Harris; puncak
dua-nada; interpolasi parabolik nada 50,3 Hz; aliasing 700 → 324 Hz;
deteksi QRS 12 detak / 72 BPM). Build dengan `gcc -std=c11 -O2
-Wall -Wextra -Wpedantic`: **0 warning**. `ctest` (CMake):
2/2 tes lolos.

## Keputusan desain

- **C11 tanpa dependensi.** Satu-satunya pustaka luar adalah libm.
  Kompleks didefinisikan sendiri (`complex.h`, `static inline`) agar
  inti FFT eksplisit dan portabel, bukan karena `<complex.h>` kurang.
- **Jendela periodik (DFT), bukan simetris.** Dengan definisi periodik,
  nada tepat di bin terbaca tepat setelah koreksi amplitudo, dan ENBW
  Hann tepat 1,50 — angka yang dikutip Harris.
- **Interpolasi parabolik pada log-magnitudo.** Mengikuti Smith
  (Spectral Audio Signal Processing): fitting parabola ke magnitudo
  dB kira-kira dua kali lebih akurat daripada magnitudo linear, dan
  hampir eksak untuk jendela Hann (galat terukur 0,016 bin untuk nada
  50,3 Hz; jendela rectangular tetap kasar, 0,166 bin).
- **Magnitudo single-sided dengan koreksi amplitudo.** `spectrum_*`
  mengembalikan amplitudo dalam satuan sinyal (nada amp A terbaca A),
  bukan |X| mentah — supaya angka yang dibaca langsung bermakna.
  Bin DC dan Nyquist tidak dikali 2 (tidak punya cerminan).
- **Detektor QRS disederhanakan dengan sadar.** Pipeline Pan–Tompkins
  (turunan → kuadrat → integrasi jendela 150 ms → ambang 25% →
  refrakter 200 ms) plus refinement ke puncak mentah ±100 ms. Ambang
  adaptif dan diskriminasi gelombang T dari paper asli tidak
  diimplementasikan — lihat Keterbatasan.
- **PRNG deterministik (xorshift64\* + Box–Muller) dengan seed
  eksplisit.** Semua sinyal berderau dapat direproduksi persis; tidak
  ada `rand()` global yang keadaan awalnya tersembunyi.
- **Format CSV teks sederhana.** Satu sampel per baris + header
  `# fs=...`. Sengaja lambat untuk file besar (I/O teks mendominasi
  waktu wall pada N = 2²⁰) demi keterbacaan dan kemudahan debug;
  format biner berada di luar scope lab ini.

## Keterbatasan & pengembangan

- **FFT hanya radix-2.** N harus pangkat dua; sinyal bukan pangkat dua
  di-zero-pad ke pangkat dua berikutnya (atau dipotong bila `--n`
  lebih kecil). DFT naif menerima N berapa pun tetapi O(N²).
- **Tidak ada zero-padding untuk interpolasi spektral.** Estimasi
  frekuensi hanya mengandalkan interpolasi parabolik tiga bin; untuk
  nada yang sangat dekat (beda < 1 bin) tetap tidak terpisahkan.
- **Detektor QRS hanya tervalidasi pada ECG sintetik** dari lab ini
  (termasuk dengan derau 0,05 mV). Belum diuji pada rekaman nyata
  (mis. MIT-BIH): ambang adaptif, diskriminasi gelombang T, dan
  penanganan baseline wander dari Pan–Tompkins (1985) belum ada.
- **ECG sintetik adalah model morfologi**, bukan simulasi fisiologis:
  tidak ada variabilitas detak (HRV), aritmia, atau artefak gerak.
- **I/O CSV teks** menjadi bottleneck pada N besar; format biner
  (mis. `.f32` raw) akan mempercepat 10–100×.
- **Belum ada STFT/spektrogram.** Analisis di sini satu frame penuh;
  sinyal non-stasioner (chirp) hanya terlihat sebagai sebaran energi,
  bukan lintasan frekuensi terhadap waktu — itu pekerjaan STFT.
- Ide pengembangan: FFT radix campuran (Bluestein untuk N sembarang),
  jendela Kaiser / Blackman–Harris, STFT + spektrogram ASCII,
  detektor QRS lengkap Pan–Tompkins dengan ambang adaptif, dan
  pembacaan file WAV.

## Lisensi

MIT — lihat [LICENSE](LICENSE).
