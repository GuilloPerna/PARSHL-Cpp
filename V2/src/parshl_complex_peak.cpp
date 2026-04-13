// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "parshl_complex_peak.hpp"
#include <cmath>
#include <algorithm>

namespace parshl {

// SAIL source: NO EQUIVALENT IN SAIL LISTING
// Paper: Serra & Smith, ICMC 1990 §5 — "Step 4: Complex Peak Interpolation"
//
// Evaluates the parabola fitted through Y[-1]=Ym1, Y[0]=Y0, Y[+1]=Yp1 at
// the fractional offset p (as returned by Qinterp on the dB spectrum).
// Used by InterpComplexAtPeak to interpolate both Re and Im FFT components
// at the true sub-bin peak location.
//
// The SAIL listing does NOT implement this step. Amplitude is taken as the raw
// dB value at the integer bin (XmagDB[MaxLoc]); phase is never extracted.
// See docs/COMPLEX_PEAK_AUDIT.md §3 for full analysis.
double ParabolaEval(double Ym1, double Y0, double Yp1, double p) noexcept {
  // Fit parabola through points at x=-1,0,+1:
  // y(x) = a x^2 + b x + c
  // c = Y0
  // a = (Ym1 + Yp1 - 2*Y0)/2
  // b = (Yp1 - Ym1)/2
  const double a = 0.5 * (Ym1 + Yp1 - 2.0*Y0);
  const double b = 0.5 * (Yp1 - Ym1);
  const double c = Y0;
  return (a*p*p + b*p + c);
}

// SAIL source: NO EQUIVALENT IN SAIL LISTING
// Paper: Serra & Smith, ICMC 1990 §5 — Step 4
//
// Evaluates the complex FFT at the interpolated peak location k+p:
//   outRe = ParabolaEval(Re[k-1], Re[k], Re[k+1], p)
//   outIm = ParabolaEval(Im[k-1], Im[k], Im[k+1], p)
//
// Input layout: SAIL !FFA interleaved convention (Parshl-source.txt L1703):
//   S[2*bin-1] = Re(bin),  S[2*bin] = Im(bin)
// Index clamping matches SAIL FindPeaks boundary guards (Parshl-source.txt L396):
//   (MaxLoc-1) MAX I1  /  (MaxLoc+1) MIN I2
//
// Called from parshl_audio_peaks.cpp for each detected partial (non-legacy mode).
// D1: interpolated magnitude via parabolic eval on Re and Im.
// D2: interpolated phase via atan2(Im(k+p), Re(k+p)).
// p comes from Qinterp (dB domain, PASO A) — not recomputed here.
void InterpComplexAtPeak(
  const sail::Array1D<double>& S,
  long Nspec,
  long k,
  double p,
  double& outRe,
  double& outIm
) {
  // Clamp indices to valid range (archaeology: match the way FindPeaks clamps its neighbors).
  long km1 = std::max<long>(1, k - 1);
  long kp1 = std::min<long>(Nspec, k + 1);

  auto Re = [&](long bin) { return S[2*bin - 1]; };
  auto Im = [&](long bin) { return S[2*bin]; };

  outRe = ParabolaEval(Re(km1), Re(k), Re(kp1), p);
  outIm = ParabolaEval(Im(km1), Im(k), Im(kp1), p);
}

} // namespace parshl
