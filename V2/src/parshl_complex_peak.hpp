// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
#include "sail_array.hpp"

namespace parshl {

// Given y(-1), y(0), y(+1) and peak offset p in [-1,1],
// return y(p) on the parabola fitted through the 3 points.
// (Standard 3-point parabolic interpolation; archaeology-friendly.)
[[nodiscard]] double ParabolaEval(double Ym1, double Y0, double Yp1, double p) noexcept;

// Interpolate complex spectrum (interleaved Re/Im) at bin k with offset p.
// S is 1-based array length 2*Nspec, with bin i stored at:
//   Re = S[2*i-1], Im = S[2*i]
void InterpComplexAtPeak(
  const sail::Array1D<double>& S,
  long Nspec,
  long k,       // integer bin index in [1..Nspec]
  double p,     // offset in [-1,1]
  double& outRe,
  double& outIm
);

} // namespace parshl
