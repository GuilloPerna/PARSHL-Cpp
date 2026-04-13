# PARSHL V2 — Archaeological Superset

**Port and extensions author:** Guillermo Perna, 2026  
**Original:** Julius O. Smith III & Xavier Serra, Stanford CCRMA, 1985  
**Criterion:** Mathematical and archaeological correctness per Smith & Serra (1987, 1990).

See `PARSHL_V2_Validation_AUDIT2_EN.pdf` for the full validation report.

---

## What V2 is

A strict superset of V1. Implements all features that Smith & Serra designed for
PARSHL but that the 1985 SAIL environment never executed:

| Feature | Flag | What it adds over V1 |
|---------|------|----------------------|
| D1/D2/D3 | active by default | Spectrally interpolated amplitude and phase at oscillator birth |
| D4 | `--reverse-analysis` | Reverse-time analysis |
| D5 | active by default | Recycling of dead oscillator slots |
| D6 | `--damax-observe` / `--damax-gate` | DAmax amplitude change monitoring |
| D7 | `--dfmax-vscale=N` | DFmax with velocity term |
| R1 | `--out-partials PREFIX` | Per-frame amp/frq export to WAV |

`--legacy` disables all of the above and reproduces V1 exactly (verified bit-identical).

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
# From this directory (V2/)
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

### Optional flags (V2)

| Flag | Description |
|------|-------------|
| `--legacy` | Reproduce V1/SAIL 1985 exactly (disables D1–D7) |
| `--reverse-analysis` | Reverse-time analysis (D4) |
| `--dfmax-vscale=N` | Scale DFmax, e.g. `--dfmax-vscale=4.0` (D7) |
| `--damax-observe=N` | Count stops where DAmax > N dB (D6) |
| `--damax-gate=N` | Active amplitude gate > N dB change (D6) |
| `--maxoscs=N` | Maximum simultaneous oscillators |
| `--out-partials PREFIX` | Export amp/frq: PREFIX.amp.wav + PREFIX.frq.wav (R1) |
| `--synth` | Resynthesize after analysis |
| `--out-wav <path>` | Save synthesis to float32 WAV |

---

## CLI — Resynthesize the canonical signals

### V2 default mode (SDR-optimal parameters)

```bash
BIN=./build/parshl_audio_peaks

# Flute A5 (~3 s, mono)  SDR=-1.427 dB
$BIN ../wav/flute-A5.wav \
     2048 256 509 -60 21.533203 0 1 0 \
     --maxoscs=32 --synth --out-wav wav_resynthesis/flute-A5_V2.wav

# Piano C4 (~36 s, stereo → channel 0)  SDR=-0.254 dB
$BIN ../wav/Piano-C4.wav \
     2048 256 6197 -60 21.533203 0 1 0 \
     --maxoscs=12 --synth --out-wav wav_resynthesis/Piano-C4_V2.wav

# Tamtam (~66 s, stereo → channel 0)  SDR=-3.291 dB
$BIN ../wav/tamtam.wav \
     2048 512 5714 -60 21.533203 0 1 0 \
     --maxoscs=16 --synth --out-wav wav_resynthesis/tamtam_V2.wav

# Female speech (~3.6 s, mono)  SDR=-2.408 dB
$BIN ../wav/female-speech.wav \
     2048 512 309 -60 21.533203 0 1 0 \
     --maxoscs=16 --synth --out-wav wav_resynthesis/female-speech_V2.wav
```

### Legacy mode (reproduces V1 exactly)

```bash
# Add --legacy to any command to reproduce V1 bit-identical
$BIN ../wav/flute-A5.wav \
     2048 256 508 -60 21.53 0 1 0 \
     --legacy --synth --out-wav flute-A5_legacy.wav
```

### Verified canonical parameters

Optimal parameters found by SDR grid search over Nfft, hop and maxoscs.

| Signal | Nfft | hop | n_frames | ThreshDB | MinSepHz | maxoscs | SDR |
|--------|------|-----|----------|----------|----------|---------|-----|
| `flute-A5.wav` | 2048 | 256 | 509 | -60 | 21.533203 | 32 | −1.427 dB |
| `Piano-C4.wav` | 2048 | 256 | 6197 | -60 | 21.533203 | 12 | −0.254 dB |
| `tamtam.wav` | 2048 | 512 | 5714 | -60 | 21.533203 | 16 | −3.291 dB |
| `female-speech.wav` | 2048 | 512 | 309 | -60 | 21.533203 | 16 | −2.408 dB |

---

## Folder contents

```
V2/
├── README.md                            ← this file
├── CMakeLists.txt                       ← build system
├── PARSHL_V2_Validation_AUDIT2_EN.pdf   ← validation report
├── src/                                 ← C++20 source code
│   ├── parshl_audio_peaks.cpp           ← main pipeline + flags D1–D7
│   ├── parshl_tracker.cpp               ← UpdateMap + D5 recycling + D6 + D7
│   ├── parshl_synthesize.cpp            ← McAulay-Quatieri synthesis
│   ├── parshl_complex_peak.cpp          ← D1–D3: InterpComplexAtPeak
│   ├── parshl_stft.cpp                  ← zero-phase STFT
│   ├── parshl_findpeaks.cpp             ← peak detection + Qinterp
│   ├── parshl_findpartials.cpp          ← FindPartials
│   ├── parshl_getwin.cpp                ← windows
│   ├── parshl_qinterp.cpp               ← parabolic sub-bin interpolation
│   ├── parshl_dbffa.cpp                 ← dBscl magnitude
│   ├── parshl_partials_writer.cpp       ← R1: oscillator export
│   ├── io_sndfile.cpp                   ← audio I/O
│   ├── fftw_wrap.cpp                    ← FFTW3 wrapper
│   ├── sail_array.hpp                   ← 1-indexed SAIL array emulation
│   └── ...
└── wav_resynthesis/
    ├── flute-A5_V2.wav      (float32, 44100 Hz, ~3 s)
    ├── Piano-C4_V2.wav      (float32, 44100 Hz, ~36 s)
    └── tamtam_V2.wav        (float32, 44100 Hz, ~66 s)
```
