# PARSHL V3 — Experimental Layer

**Port and extensions author:** Guillermo Perna, 2026  
**Original:** Julius O. Smith III & Xavier Serra, Stanford CCRMA, 1985  
**Criterion:** Verified improvement in tracking and synthesis metrics over the V2 baseline.

See `PARSHL_V3_Experimental_AUDIT2_EN.pdf` for the full report.

---

## What V3 is

A strict superset of V2. Without flags, V3 produces bit-identical results to V2.

The **recommended production combination** was found by grid search
(108 combinations × 3 signals). Verified improvements over V2:

- −24% tracking stops on Piano-C4
- +41% avg_lifetime on Piano-C4
- +4.65 dB SDR on flute-A5
- Elimination of the SAIL amplitude bias

```bash
# Recommended production flags
--sigscl-analytic --oscphs-double --closest-osc-free
```

---

## Dependencies

### Install on Ubuntu/Debian

```bash
sudo apt update
sudo apt install -y cmake build-essential libfftw3-dev libsndfile1-dev pkg-config
```

| Package | Min version | Purpose |
|---------|-------------|---------|
| `cmake` | 3.20 | Build system |
| `build-essential` | — | GCC/G++ with C++20 support (GCC >= 10) |
| `libfftw3-dev` | 3.x | Fast FFT |
| `libsndfile1-dev` | 1.x | WAV audio read/write |
| `pkg-config` | — | Library detection for CMake |

### Verify before building

```bash
cmake --version        # >= 3.20
gcc --version          # >= 10
pkg-config --modversion sndfile   # e.g. 1.2.2
pkg-config --modversion fftw3     # e.g. 3.3.10
```

---

## Build

```bash
# From this directory (V3/)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The binary is at `build/parshl_audio_peaks`.

---

## Usage

```
./build/parshl_audio_peaks <file.wav> <Nfft> <hop> <n_frames> <ThreshDB> \
    <MinSepHz> <Hyst> <MinWid> <tracePartials> \
    [flags] [--synth] [--out-wav <output.wav>]
```

### Positional parameters (in order)

| # | Parameter | Description |
|---|-----------|-------------|
| 1 | `file.wav` | Input (mono or stereo; stereo uses channel 0) |
| 2 | `Nfft` | FFT size (2048 or 4096) |
| 3 | `hop` | Frame step in samples (256 or 128) |
| 4 | `n_frames` | Frames to analyse (see canonical table) |
| 5 | `ThreshDB` | Detection threshold in dBscl (e.g. -60) |
| 6 | `MinSepHz` | Minimum peak separation in Hz (see canonical table) |
| 7 | `Hyst` | Tracker hysteresis. Use **0** |
| 8 | `MinWid` | Minimum peak width. Use **1** |
| 9 | `tracePartials` | 1 = verbose, 0 = silent |

### Production E-feature flags (recommended)

| Flag | Experiment | Measured effect |
|------|------------|----------------|
| `--sigscl-analytic` | E1 | Analytic OLA normalisation 1/Σw[n]. Eliminates SAIL amplitude bias. |
| `--oscphs-double` | E4 | OscPhs in double precision. Eliminates phase drift. +4.65 dB SDR flute. |
| `--closest-osc-free` | E2 | Recycles the oscillator with the closest frequency. Fewer assignment errors. |

### All available flags

| Flag | Description |
|------|-------------|
| `--legacy` | Reproduce V1/SAIL 1985 exactly (disables everything) |
| `--sigscl-analytic` | E1: analytic SigScl |
| `--oscphs-double` | E4: double-precision OscPhs |
| `--closest-osc-free` | E2: closest-frequency oscillator recycling |
| `--dfmax-vscale=N` | E2/D7: DFmax multiplier |
| `--reverse-analysis` | D4: reverse-time analysis |
| `--damax-observe=N` | D6: observe DAmax stops |
| `--damax-gate=N` | D6: active amplitude gate |
| `--maxoscs=N` | Maximum simultaneous oscillators |
| `--out-partials PREFIX` | R1: export per-frame amp/frq to WAV |
| `--synth` | Resynthesize after analysis |
| `--out-wav <path>` | Save synthesis to float32 WAV |

---

## CLI — Resynthesize the canonical signals

### Production mode (SDR-optimal parameters)

```bash
BIN=./build/parshl_audio_peaks
OPT="--sigscl-analytic --oscphs-double --closest-osc-free"

# Flute A5 (~3 s, mono)  SDR=+0.650 dB
$BIN ../wav/flute-A5.wav \
     2048 256 509 -60 21.533203 0 1 0 \
     --maxoscs=10 $OPT --synth --out-wav wav_resynthesis/flute-A5_V3.wav

# Piano C4 (~36 s, stereo → channel 0)  SDR=+1.434 dB
$BIN ../wav/Piano-C4.wav \
     2048 256 6197 -60 21.533203 0 1 0 \
     --maxoscs=16 $OPT --synth --out-wav wav_resynthesis/Piano-C4_V3.wav

# Tamtam (~66 s, stereo → channel 0)  SDR=-1.076 dB
$BIN ../wav/tamtam.wav \
     2048 512 5714 -60 21.533203 0 1 0 \
     --maxoscs=32 $OPT --synth --out-wav wav_resynthesis/tamtam_V3.wav

# Female speech (~3.6 s, mono)  SDR=-0.333 dB
$BIN ../wav/female-speech.wav \
     2048 512 309 -60 21.533203 0 1 0 \
     --maxoscs=24 $OPT --synth --out-wav wav_resynthesis/female-speech_V3.wav
```

### V2 baseline mode (no production flags)

```bash
# No flags = bit-identical to V2 result
$BIN ../wav/flute-A5.wav \
     2048 256 508 -60 21.53 0 1 0 \
     --synth --out-wav flute-A5_V2baseline.wav
```

### Legacy mode (identical to V1)

```bash
$BIN ../wav/flute-A5.wav \
     2048 256 508 -60 21.53 0 1 0 \
     --legacy --synth --out-wav flute-A5_V1legacy.wav
```

### Verified canonical parameters

Optimal parameters found by SDR grid search over Nfft, hop and maxoscs.

| Signal | Nfft | hop | n_frames | ThreshDB | MinSepHz | maxoscs | SDR |
|--------|------|-----|----------|----------|----------|---------|-----|
| `flute-A5.wav` | 2048 | 256 | 509 | -60 | 21.533203 | 10 | +0.650 dB |
| `Piano-C4.wav` | 2048 | 256 | 6197 | -60 | 21.533203 | 16 | +1.434 dB |
| `tamtam.wav` | 2048 | 512 | 5714 | -60 | 21.533203 | 32 | −1.076 dB |
| `female-speech.wav` | 2048 | 512 | 309 | -60 | 21.533203 | 24 | −0.333 dB |

---

## Negative experiments (do not use)

| Flag | Reason for closure |
|------|--------------------|
| `--dbspec-refine-amp` | Breaks the (freq, amp) invariant from the same interpolation point p → corrupts OLA |
| `--dbspec-match` | Only viable at factor=0.20 and even then costs −0.43 dB SDR on Piano-C4 |
| `--grs-groups=N` | Band centroids ≠ spectral peaks → incoherent grouping, erroneous partial assignments |

---

## Folder contents

```
V3/
├── README.md                              ← this file
├── CMakeLists.txt                         ← build system
├── PARSHL_V3_Experimental_AUDIT2_EN.pdf   ← experimental report
├── src/                                   ← C++20 source code
│   ├── parshl_audio_peaks.cpp             ← main pipeline + D1–D7 + E1–E6
│   ├── parshl_tracker.cpp                 ← tracker + D5/D6/D7/E2
│   ├── parshl_synthesize.cpp              ← synthesis + E4 (oscphs-double)
│   ├── parshl_complex_peak.cpp            ← D1–D3: InterpComplexAtPeak
│   ├── parshl_stft.cpp                    ← zero-phase STFT
│   ├── parshl_findpeaks.cpp               ← peak detection
│   ├── parshl_findpartials.cpp            ← FindPartials + E1
│   ├── parshl_getwin.cpp                  ← windows
│   ├── parshl_qinterp.cpp                 ← parabolic sub-bin interpolation
│   ├── parshl_dbffa.cpp                   ← dBscl magnitude
│   ├── parshl_partials_writer.cpp         ← R1: oscillator export
│   ├── io_sndfile.cpp                     ← audio I/O
│   ├── fftw_wrap.cpp                      ← FFTW3 wrapper
│   ├── sail_array.hpp                     ← 1-indexed SAIL array emulation
│   └── ...
└── wav_resynthesis/
    ├── flute-A5_V3.wav      (float32, 44100 Hz, ~3 s)
    ├── Piano-C4_V3.wav      (float32, 44100 Hz, ~36 s)
    └── tamtam_V3.wav        (float32, 44100 Hz, ~66 s)
```
