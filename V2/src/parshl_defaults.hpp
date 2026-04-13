// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once

// Canonical PARSHL defaults reconstructed from original SAIL source.
// DO NOT CHANGE unless justified by SAIL listing evidence.

namespace parshl_defaults {

constexpr int   Nfft  = 1024;
constexpr int   Nx    = 676;
constexpr int   Nhop  = Nx / 2;
constexpr int   Nh    = 0;
constexpr int   Ndec  = 1;

enum class WindowType {
Hamming,
Hann,
Rectangular,
Triangular
};

constexpr WindowType WinType = WindowType::Hamming;

constexpr bool DoSynth      = true;
constexpr bool DoFlt        = false;
constexpr bool SwapOut      = false;
constexpr bool InstantRise  = false;

constexpr int Trace = 0;

} // namespace parshl_defaults

