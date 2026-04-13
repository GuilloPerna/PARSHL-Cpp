// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port, archaeological extensions (V2), and experimental
//   improvements (V3): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// src/parshl_findpeaks.hpp
#pragma once
#include "sail_array.hpp"
#include "parshl_globals.hpp"
#include "parshl_qinterp.hpp"
#include <cstdint>

namespace parshl {

// ORIG (SAIL):
// INTEGER PROCEDURE FindPeaks(REAL ARRAY X,Peaks,PeakLocs;
//   REAL Thresh(RealNoSet),Hyst(RealNoSet);
//   INTEGER MaxPeaks(IntNoSet),MinWidth(IntNoSet),BlankWidth(IntNoSet),
//           I1(IntNoSet),I2(IntNoSet));
//
// Notes (SAIL):
// - X is clobbered (local maxima are removed by "clobbering" to a valley value).
// - Peaks[]: peak height (MaxVal) in same units as X (here, dB).
// - PeakLocs[]: MaxLoc + Qinterp(...) where Qinterp returns offset in [-1,1] (clipped).
[[nodiscard]] int FindPeaks(
  sail::Array1D<double>& X,
  sail::Array1D<double>& Peaks,
  sail::Array1D<double>& PeakLocs,
  double Thresh = RealNoSet,
  double Hyst = RealNoSet,
  std::int64_t MaxPeaks = IntNoSet,
  std::int64_t MinWidth = IntNoSet,
  std::int64_t BlankWidth = IntNoSet,
  std::int64_t I1 = IntNoSet,
  std::int64_t I2 = IntNoSet
);

} // namespace parshl
