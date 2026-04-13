// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
#include <cstdint>
#include "sail_array.hpp"          // tu Array1D
#include "parshl_findpeaks.hpp"    // FindPeaks
#include "parshl_qisort.hpp"       // QIsort

namespace parshl {

// FindPartials: archaeological wrapper around FindPeaks.
// Inputs:
//  - Xdb: spectrum in dB (Array1D, 1..Nspec)
//  - Fs, Nfft, Fc1, Fc2, MinSep: PARSHL parameters
//  - Thresh, Hyst, MinWid: FindPeaks parameters
// Outputs:
//  - LinAmpDB: magnitudes in dB per line (1..NpartialsFound)
//  - LinFrq: frequencies in Hz per line (1..NpartialsFound)
// Returns: NpartialsFound
[[nodiscard]] std::int64_t FindPartials(
  const sail::Array1D<double>& Xdb,
  sail::Array1D<double>& LinAmpDB,
  sail::Array1D<double>& LinFrq,
  double Fs,
  std::int64_t Nfft,
  double Fc1,
  double Fc2,
  double MinSep,
  double Thresh,
  double Hyst,
  std::int64_t NpartialsReq,
  std::int64_t MinWid,
  bool trace=false
);

} // namespace parshl
