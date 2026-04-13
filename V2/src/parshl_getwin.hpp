// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
#include "sail_array.hpp"
#include <cstdint>
#include <stdexcept>

namespace parshl {

// WinType enum — V1 correction: reordered to match SASP/JOS convention.
// SAIL original had Hamming=3, Hanning=5 (swapped vs standard usage).
// Corrected order: Hanning=3, Hamming=4, GenHamming=5.
enum WinType : int {
  Rectangular = 1,
  Triangular  = 2,
  Hanning     = 3,
  Hamming     = 4,
  GenHamming  = 5,
  Kaiser      = 6,
  Chebyshev   = 7
};

// Port of GetWin(W, Wtype, Nw, P3, P4).
// - W is 1-based [1:Nw]
// - P3/P4 meaning matches the SAIL prompts: Kaiser uses stopband rejection dB (P3), GenHamming uses alpha (P3),
///  Chebyshev uses stopband rejection dB (P3) and/or transition width as negative P4 trick in SAIL. :contentReference[oaicite:3]{index=3}
void GetWin(sail::Array1D<double>& W, int Wtype, int Nw, double P3 = -1.0, double P4 = -1.0);

} // namespace parshl
