# PARSHL V1 — Faithful Port of SAIL 1985

**Port author:** Guillermo Perna, 2026  
**Original:** Julius O. Smith III & Xavier Serra, Stanford CCRMA, 1985  
**Criterion:** Structural, algorithmic and numerical fidelity to the original SAIL.

See `PARSHL_V1_Validation_AUDIT2_EN.pdf` for the full validation report.

---

## What V1 is

A 1:1 port of the original PARSHL system written in SAIL/ALGOL for the DEC-10.
No modernisations or algorithmic improvements. Five platform corrections were applied
to adapt the code to the C++20/Linux environment (see PDF report).

The `parshl_audio_peaks` binary analyses a WAV file, detects sinusoidal partials
frame by frame, and resynthesises using the McAulay-Quatieri oscillator bank method.

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
| `libfftw3-dev` | 3.x | Fast FFT (replaces SAIL !FFA/!FFS) |
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
# From this directory (V1/)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The binary is at `build/parshl_audio_peaks`.  
Build takes 10–30 seconds.

### Verify the build

```bash
ls -lh build/parshl_audio_peaks
```

---

## Usage

```
./build/parshl_audio_peaks <file.wav> <Nfft> <hop> <n_frames> <ThreshDB> \
    <MinSepHz> <Hyst> <MinWid> <tracePartials> \
    [--legacy] [--synth] [--out-wav <output.wav>]
```

### Positional parameters (in order)

| # | Parameter | Type | Description |
|---|-----------|------|-------------|
| 1 | `file.wav` | path | Input signal (mono or stereo; stereo uses channel 0) |
| 2 | `Nfft` | integer | FFT size. Recommended: 2048 (short signals) or 4096 (Piano) |
| 3 | `hop` | integer | Frame step in samples. Recommended: 256 or 128 |
| 4 | `n_frames` | integer | Number of frames to analyse (see canonical table below) |
| 5 | `ThreshDB` | float | Peak detection threshold in dBscl. Recommended: -60 or -70 |
| 6 | `MinSepHz` | float | Minimum separation between peaks in Hz (see canonical table) |
| 7 | `Hyst` | integer | Tracker hysteresis. Use **0** |
| 8 | `MinWid` | integer | Minimum peak width in bins. Use **1** |
| 9 | `tracePartials` | 0/1 | 1 = print partials per frame (verbose), 0 = silent |

### Optional flags

| Flag | Description |
|------|-------------|
| `--legacy` | SAIL 1985 faithful mode. Recommended to document V1 intent. |
| `--synth` | Resynthesize after analysis. Without this flag, analysis only. |
| `--out-wav <path>` | Save synthesis to a float32 WAV file. |
| `--maxoscs=N` | Maximum simultaneous oscillators |
| `--normalize-peak=<dBFS>` | N1: post-synthesis peak normalization. Example: `--normalize-peak=-1.0`. Not SAIL-faithful. |

---

## CLI — Resynthesize the canonical signals

Run from the `V1/` directory:

```bash
BIN=./build/parshl_audio_peaks

# Flute A5 (~3 s, mono)  SDR=-2.162 dB
$BIN ../wav/flute-A5.wav \
     1024 256 513 -60 43.066406 0 1 0 676 \
     --maxoscs=12 --synth --out-wav wav_resynthesis/flute-A5_V1.wav

# Piano C4 (~36 s, stereo → channel 0)  SDR=-2.520 dB
$BIN ../wav/Piano-C4.wav \
     1024 256 6201 -60 43.066406 0 1 0 1024 \
     --maxoscs=16 --synth --out-wav wav_resynthesis/Piano-C4_V1.wav

# Tamtam (~66 s, stereo → channel 0)  SDR=-2.009 dB
$BIN ../wav/tamtam.wav \
     1024 256 11432 -60 43.066406 0 1 0 676 \
     --maxoscs=16 --synth --out-wav wav_resynthesis/tamtam_V1.wav

# Female speech (~3.6 s, mono)  SDR=-2.414 dB
$BIN ../wav/female-speech.wav \
     1024 256 622 -60 43.066406 0 1 0 1024 \
     --maxoscs=16 --synth --out-wav wav_resynthesis/female-speech_V1.wav
```

### Normalized WAV (convenience, non-SAIL-faithful)

```bash
# Flute A5 at −1 dBFS — prevents DAC clipping; not part of default SAIL behaviour
$BIN ../wav/flute-A5.wav \
     1024 256 513 -60 43.066406 0 1 0 676 \
     --maxoscs=12 --normalize-peak=-1.0 --synth --out-wav wav_resynthesis/flute-A5_normalized_V1.wav
```

Output WAVs are float32, 44100 Hz, mono.

### Verified canonical parameters

Optimal parameters found by SDR grid search over Nfft, hop and maxoscs.
The 10th positional argument `Nx` sets the SigScl heuristic: SigScl = 4/Nx.

| Signal | Nfft | hop | n_frames | ThreshDB | MinSepHz | Nx | maxoscs | SDR |
|--------|------|-----|----------|----------|----------|----|---------|-----|
| `flute-A5.wav` | 1024 | 256 | 513 | -60 | 43.066406 | 676 | 12 | −2.162 dB |
| `Piano-C4.wav` | 1024 | 256 | 6201 | -60 | 43.066406 | 1024 | 16 | −2.520 dB |
| `tamtam.wav` | 1024 | 256 | 11432 | -60 | 43.066406 | 676 | 16 | −2.009 dB |
| `female-speech.wav` | 1024 | 256 | 622 | -60 | 43.066406 | 1024 | 16 | −2.414 dB |

**Important:** The exact values of n_frames and MinSepHz are critical.
Use the values in this table, not approximations.

---

## Folder contents

```
V1/
├── README.md                            ← this file
├── CMakeLists.txt                       ← build system
├── PARSHL_V1_Validation_AUDIT2_EN.pdf   ← validation report
├── src/                                 ← C++20 source code
│   ├── parshl_audio_peaks.cpp           ← main pipeline
│   ├── parshl_tracker.cpp               ← partial tracking (UpdateMap)
│   ├── parshl_synthesize.cpp            ← McAulay-Quatieri synthesis
│   ├── parshl_stft.cpp                  ← STFT with zero-phase layout (C1)
│   ├── parshl_findpeaks.cpp             ← peak detection in dB
│   ├── parshl_findpartials.cpp          ← FindPartials → LinFrq/LinAmpDB
│   ├── parshl_complex_peak.cpp          ← sub-bin Re/Im interpolation
│   ├── parshl_getwin.cpp                ← windows (Hamming, Hanning, etc.)
│   ├── parshl_qinterp.cpp               ← parabolic sub-bin interpolation
│   ├── parshl_dbffa.cpp                 ← dBscl magnitude
│   ├── parshl_partials_writer.cpp       ← partial export
│   ├── io_sndfile.cpp                   ← audio I/O with ×32768 scale (C2)
│   ├── fftw_wrap.cpp                    ← FFTW3 wrapper
│   ├── sail_array.hpp                   ← 1-indexed SAIL array emulation
│   └── ...
└── wav_resynthesis/
    ├── flute-A5_V1.wav      (float32, 44100 Hz, ~3 s)
    ├── Piano-C4_V1.wav      (float32, 44100 Hz, ~36 s)
    └── tamtam_V1.wav        (float32, 44100 Hz, ~66 s)
```
