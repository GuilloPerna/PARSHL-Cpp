// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port, archaeological extensions (V2), and experimental
//   improvements (V3): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// src/parshl_qinterp.cpp
#include "parshl_qinterp.hpp"
#include "parshl_globals.hpp"
#include <cmath>
#include <iostream>

namespace parshl {

// SAIL source: INTERNAL SIMPLE REAL PROCEDURE Qinterp
// Parshl-source.txt L315–337
//
// Fits a parabola Y(X) = A*X^2 + B through three equally spaced dB values
//   Y[-1]=Ym1,  Y[0]=Y0,  Y[+1]=Yp1
// and returns the X-coordinate (offset) where the extremum is attained.
//
// Used in FindPeaks (Parshl-source.txt L396) on the dB magnitude spectrum
// (XmagDB) to compute the sub-bin frequency offset of each spectral peak.
// Step 3 of the analysis algorithm (paper §5).
//
// Inputs:
//   Ym1, Y0, Yp1 — three consecutive dB values around the peak bin
//   InteriorX    — if true, clip result X to [-1, 1] (SAIL default: TRUE)
// Returns: X ∈ [-1,1] — fractional bin offset of the parabola extremum
//
// ORIG (SAIL L328–334):
//   X <- (Yp1 - Ym1) / (2*(2*Y0 - Yp1 - Ym1))
//   IF InteriorX AND (ABS(X) > 1): X <- (IF X>0 THEN 1 ELSE -1)
//   IF Debug3 THEN PRINT(...)
double Qinterp(double Ym1, double Y0, double Yp1, bool InteriorX) {
  double denom = 2.0 * (2.0 * Y0 - Yp1 - Ym1);
  double X = (Yp1 - Ym1) / denom;

  if (InteriorX && (std::abs(X) > 1.0)) {
    std::cout << " Qinterp: Clipping analytic extremum to ";
    X = (X > 0.0) ? 1.0 : -1.0;
    std::cout << X << "\n";
  }

  if (Debug3) {
    std::cout << "Qinterp: Given Y's " << Ym1 << " " << Y0 << " " << Yp1
              << ", extremum is at X = " << X << "\n";
  }
  return X;
}

} // namespace parshl

//Note: The formula matches the SAIL listing (parabola fitted over three equally spaced points)
