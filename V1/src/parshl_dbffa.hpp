// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
#include "sail_array.hpp"

namespace parshl {

// ORIG: PROCEDURE DpyFFA(REAL ARRAY S; INTEGER N; REAL DbMin,DbMax; INTEGER PwrOfTen(0); BOOLEAN InDecibels(TRUE));
// We only port the "compute dB array" core, not plotting.
// Input S is complex interleaved (Re/Im pairs), 1-based indexing.
// Output Xdb is real, 1-based indexing, length N/2 (positive-frequency bins count depends on caller).
void DbFromInterleavedComplex(
  const sail::Array1D<double>& S,
  long N,                 // total length of S in re/im scalars (even)
  sail::Array1D<double>& Xdb
);

} // namespace parshl
