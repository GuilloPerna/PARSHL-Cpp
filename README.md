# PARSHL — C++20 Port (V1/V2/V3)

**Port and extensions author:** Guillermo Perna, 2026  
**Original SAIL program:** Julius O. Smith III, Stanford CCRMA, 1985  
**1987 paper (V2 features):** Julius O. Smith III & Xavier Serra, ICMC-87 / CMJ 14(4), 1990

---

## Package contents

```
PARSHL_Cpp/
├── README.md                                    ← this file
├── PARSHL_V1_V2_V3_Comparative_AUDIT2_EN.pdf   ← three-version comparative report
├── wav/   ← canonical test signals (shared by all versions)
├── V1/    ← Faithful port of Smith's SAIL 1985 (wav_resynthesis/ inside)
├── V2/    ← Smith & Serra 1987 paper features (wav_resynthesis/ inside)
└── V3/    ← Author's experimental layer (wav_resynthesis/ inside)
```

Each folder contains: complete source code, CMakeLists.txt, test WAVs, README with build and CLI instructions, and validation report PDF.

---

## The three versions

| Version | What it is | When to use it |
|---------|------------|----------------|
| **V1** | 1:1 port of the SAIL 1985 original. Reproduces the historical behaviour exactly. | Archaeological reference. For studying Smith's original SAIL system. |
| **V2** | V1 + features of the paper not implemented in 1985: complex spectral interpolation + phase init (D1–D3), reverse-time analysis (D4), oscillator recycling (D5), DAmax monitoring (D6), adaptive DFmax (D7). `--legacy` = V1 exact. | To hear PARSHL as designed in the 1987 paper. |
| **V3** | V2 + author's own experiments: analytic OLA normalisation (E1), closest-frequency oscillator assignment (E2), double-precision OscPhs (E4). No flags = V2 exact. | For the best available synthesis quality. |

**Inheritance:** V1 ⊂ V2 ⊂ V3 (strict supersets, bit-identical backwards)

---

## Dependencies — installation

### Supported operating system

Ubuntu 22.04 / 24.04 (or any modern Debian-based distribution).  
Verified on Ubuntu 24.04 with GCC 13.3, CMake 3.28, FFTW3 3.3.10, libsndfile 1.2.2.

### Required packages

```bash
sudo apt update
sudo apt install -y \
    cmake \
    build-essential \
    libfftw3-dev \
    libsndfile1-dev \
    pkg-config
```

| Package | Min version | Purpose |
|---------|-------------|---------|
| `cmake` | 3.20 | Build system |
| `build-essential` | — | GCC/G++ C++20 compiler |
| `libfftw3-dev` | 3.x | Fast Fourier Transform (replaces SAIL !FFA/!FFS) |
| `libsndfile1-dev` | 1.x | WAV audio I/O |
| `pkg-config` | — | Library detection for CMake |

### Verify installation

```bash
cmake --version          # must be 3.20 or higher
gcc --version            # must be GCC 10 or higher
pkg-config --modversion sndfile   # e.g. 1.2.2
pkg-config --modversion fftw3     # e.g. 3.3.10
```

---

## Build (same for V1, V2 and V3)

```bash
# From the chosen version folder, e.g. V3:
cd V3
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The binary is at `build/parshl_audio_peaks`.  
Build takes 10–30 seconds depending on hardware.

---

## Comparative report

`PARSHL_V1_V2_V3_Comparative_AUDIT2_EN.pdf` compares all three versions:
- Feature table per version
- True tracking baseline (full signals, 2026-04-05)
- V3 grid search results
- JOS tests T1/T2/T3

---

## Canonical test signals

The four signals are in `wav/` (shared root directory; identical for all three versions).

| Signal | Duration | Format |
|--------|----------|--------|
| `flute-A5.wav` | ~3 s | mono, 44100 Hz, 16-bit PCM |
| `Piano-C4.wav` | ~36 s | stereo, 44100 Hz, 16-bit PCM |
| `tamtam.wav` | ~66 s | stereo, 44100 Hz, 16-bit PCM |
| `female-speech.wav` | ~3.6 s | mono, 44100 Hz, 16-bit PCM |

**Note on stereo files:** `Piano-C4.wav` and `tamtam.wav` are stereo.
PARSHL automatically analyses channel 0 (left).

---

## Resynthesis audio samples

Each version folder contains a `wav_resynthesis/` directory with
additive synthesis output for all four signals, generated with
SDR-optimal parameters:

| Signal | V1 SDR | V2 SDR | V3 SDR |
|--------|--------|--------|--------|
| `flute-A5`      | −2.16 dB | −1.43 dB | −0.704 dB |
| `Piano-C4`      | −2.52 dB | −0.25 dB | **+0.484 dB** |
| `tamtam`        | −2.01 dB | −3.29 dB | −2.955 dB |
| `female-speech` | −2.414 dB | −2.408 dB | −1.978 dB |

V1 uses Nfft=1024 (historical SAIL parameters, Smith 1985).  
V2 uses Nfft=2048 with complex spectral interpolation (Smith & Serra 1987,
no phase corrections).  
V3 uses Nfft=2048 + analytic DFT normalisation (`--sigscl-analytic`, E1 — corrected
to `SigScl = 2/Σw`) + double-precision OscPhs (`--oscphs-double`, E4) +
closest-frequency assignment (`--closest-osc-free`, E2).

**Amplitude notes:**
- **PARSHL original has no output normalization or overshoot protection.**
  The original SAIL source (`Parshl-source.txt`) initializes `Maxamp ← 0` (L942) and
  never updates it from the synthesis buffer — it is passed to `WriteH()` (L1795) only
  as file-header metadata. The `Sando()` output routine writes samples raw, with no
  clamping. **By default, this port does not prevent overshoot either: synthesized
  samples may exceed ±1.0 and will clip on playback or when saved as integer PCM.**
  V1 and V2 reproduce this SAIL behaviour exactly by default.
- **V1/V2:** the SAIL heuristic `SigScl = 4/Nfft` is accurate for Hann windows but
  introduces a **+0.67 dB** positive bias for Hamming windows. This causes overshoot
  on `flute-A5` (+2.47 dBFS / +2.77 dBFS for V1/V2; 2303 / 4129 samples ≥ 1.0).
  Normalized alternatives at −1 dBFS are provided as `flute-A5_normalized_V1.wav` /
  `flute-A5_normalized_V2.wav`. The optional flag `--normalize-peak=-1.0` prevents
  overshoot on any signal, but is **not** SAIL-faithful; it is absent from the original PARSHL.
- **V3:** `--sigscl-analytic` applies `SigScl = 2/Σw` which is window-independent and
  produces RMS within **−1 dB** of the original across all four test signals (no overshoot).

The negative SDR on `tamtam` across all versions reflects the theoretical
limit of sinusoidal additive synthesis for inharmonic percussive signals:
the broad-band noise component of the attack cannot be represented by a
finite sum of sinusoids, regardless of parameter settings.
