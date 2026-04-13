// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// src/parshl_findpeaks.cpp
#include "parshl_findpeaks.hpp"

#include <algorithm>
#include <cstdint>

namespace parshl {

// --- local helpers (SAIL ArrMin / ArrMax analogs) ---
[[nodiscard]] static inline std::int64_t clamp_i64(std::int64_t v, std::int64_t lo, std::int64_t hi) noexcept {
  return std::max(lo, std::min(hi, v));
}

static std::int64_t ArrMin(double& outMin, const sail::Array1D<double>& X,
                           std::int64_t I1, std::int64_t I2) {
  // Returns index of min (not used by SAIL here, but convenient).
  std::int64_t idx = I1;
  double mn = X[I1];
  for (std::int64_t i = I1 + 1; i <= I2; ++i) {
    if (X[i] < mn) { mn = X[i]; idx = i; }
  }
  outMin = mn;
  return idx;
}

static std::int64_t ArrMax(double& outMax, const sail::Array1D<double>& X,
                           std::int64_t I1, std::int64_t I2) {
  // Returns index of max, sets outMax.
  std::int64_t idx = I1;
  double mx = X[I1];
  for (std::int64_t i = I1 + 1; i <= I2; ++i) {
    if (X[i] > mx) { mx = X[i]; idx = i; }
  }
  outMax = mx;
  return idx;
}

// SAIL source: INTEGER PROCEDURE FindPeaks
// Parshl-source.txt L344–L427
//
// Iteratively locates amplitude peaks in the dB magnitude spectrum X[I1:I2].
// Algorithm per peak:
//   1. Find global maximum via ArrMax in the remaining search range      (L387–L388)
//   2. Compute sub-bin location via Qinterp on three neighbouring dB values (L396)
//   3. Record raw dB amplitude at the integer bin                          (L397)
//   4. Expand blank interval [M1,M2] left/right while spectrum is flat      (L404–L414)
//      (Hyst threshold prevents re-detection of the same local maximum)
//   5. Set X[M1:M2] = ClobVal (minimum of the flanking region)             (L413)
//   6. Reject peak if too narrow or at I1/I2 boundary ("throw it back")    (L416–L423)
//   7. Repeat until MaxPeaks found or no candidate exceeds Thresh
//
// X is modified in-place (blanked). Caller (FindPartials) saves & restores it.
// Outputs PeakLocs[] in fractional bin units; amplitude Peaks[] in dB.
int FindPeaks(
  sail::Array1D<double>& X,
  sail::Array1D<double>& Peaks,
  sail::Array1D<double>& PeakLocs,
  double Thresh,
  double Hyst,
  std::int64_t MaxPeaks,
  std::int64_t MinWidth,
  std::int64_t BlankWidth,
  std::int64_t I1,
  std::int64_t I2
) {
  // Bounds
  const std::int64_t Lb = X.lo();
  const std::int64_t Ub = X.hi();

  // SAIL defaults:
  // IF I2=IntNoSet THEN I2 <- ARRINFO(X,2);
  // IF I1=IntNoSet THEN I1 <- ARRINFO(X,1);
  if (I2 == IntNoSet) I2 = Ub;
  if (I1 == IntNoSet) I1 = Lb;

  // I1 <- (Lb MAX I1 MIN Ub);  I2 <- (I1 MAX I2 MIN Ub);
  I1 = clamp_i64(I1, Lb, Ub);
  I2 = clamp_i64(I2, I1, Ub);

  // Poff <- ARRINFO(Peaks,1)-1; PLoff <- ARRINFO(PeakLocs,1)-1;
  const std::int64_t Poff  = Peaks.lo()    - 1;
  const std::int64_t PLoff = PeakLocs.lo() - 1;

  // Npeaks <- (IF MaxPeaks NEQ IntNoSet THEN MaxPeaks ELSE ARRINFO(Peaks,2)-Poff);
  std::int64_t Npeaks = (MaxPeaks != IntNoSet) ? MaxPeaks : (Peaks.hi() - Poff);

  // Npeaks <- Npeaks MIN (ARRINFO(Peaks,2)-Poff) MIN (ARRINFO(PeakLocs,2)-PLoff);
  Npeaks = std::min<std::int64_t>(Npeaks, (Peaks.hi()    - Poff));
  Npeaks = std::min<std::int64_t>(Npeaks, (PeakLocs.hi() - PLoff));

  // NdxReach <- (IF (BlankWidth NEQ IntNoSet) THEN (BlankWidth-1)/2 MAX 1 ELSE 0);
  std::int64_t NdxReach = 0;
  if (BlankWidth != IntNoSet) {
    NdxReach = (BlankWidth - 1) / 2;
    if (NdxReach < 1) NdxReach = 1;
  }

  // ArrMin/ArrMax over [I1:I2]
  double Xmin = 0.0, Xmax = 0.0;
  ArrMin(Xmin, X, I1, I2);
  ArrMax(Xmax, X, I1, I2);

  // IF Thresh=RealNoSet THEN Thresh <- Xmin;
  // IF Hyst=RealNoSet THEN Hyst <- (Xmax-Xmin)/100;
  if (Thresh == RealNoSet) Thresh = Xmin;
  if (Hyst == RealNoSet)   Hyst   = (Xmax - Xmin) / 100.0;

  // M1 <- M2 <- 0; COMMENT [m1,m2] = index interval of last peak;
  std::int64_t M1 = 0, M2 = 0;
  std::int64_t Nfound = 0;

  for (std::int64_t Peak = 1; Peak <= Npeaks; ++Peak) {
    // IF M1=I1 AND M2=I2 THEN DONE "fp";
    if (M1 == I1 && M2 == I2) break;

    // MaxLoc <- ArrMax(MaxVal,X,I1,I2);
    double MaxVal = 0.0;
    const std::int64_t MaxLoc = ArrMax(MaxVal, X, I1, I2);

    // IF MaxVal<Thresh THEN DONE "fp";
    if (MaxVal < Thresh) break;

    // Nfound <- Nfound + 1;
    Nfound++;

    // SAIL L396: PeakLocs[n+PLoff] <- MaxLoc + Qinterp(X[(MaxLoc-1) MAX I1], MaxVal, X[(MaxLoc+1) MIN I2])
    // Parabolic interpolation on dB spectrum gives sub-bin frequency offset.
    const std::int64_t L = std::max<std::int64_t>(MaxLoc - 1, I1);
    const std::int64_t R = std::min<std::int64_t>(MaxLoc + 1, I2);
    PeakLocs[Nfound + PLoff] = static_cast<double>(MaxLoc) + Qinterp(X[L], MaxVal, X[R], true);

    // SAIL L397: Peaks[Nfound+Poff] <- MaxVal   (raw dB at integer bin)
    Peaks[Nfound + Poff] = MaxVal;

    // SAIL L400: M1 <- (MaxLoc-NdxReach) MAX I1;  M2 <- (MaxLoc+NdxReach) MIN I2
    // COMMENT Now slice off peak so we don't find it again;
    // M1 <- (MaxLoc-NdxReach) MAX I1;  M2 <- (MaxLoc+NdxReach) MIN I2;
    M1 = std::max<std::int64_t>(MaxLoc - NdxReach, I1);
    M2 = std::min<std::int64_t>(MaxLoc + NdxReach, I2);

    // ArrMin(ClobVal,X,M1,M2);
    double ClobVal = 0.0;
    ArrMin(ClobVal, X, M1, M2);

    // SAIL L404–L406: expand M1 left while X[M1-1] is within Hyst of the trough
    // TmpR <- X[M1];
    // WHILE M1>I1 AND TmpR+Hyst GEQ X[M1-1] DO BEGIN TmpR <- TmpR MIN X[M1-1]; M1 <- M1-1 END;
    double TmpR = X[M1];
    while (M1 > I1 && (TmpR + Hyst) >= X[M1 - 1]) {
      TmpR = std::min(TmpR, X[M1 - 1]);
      M1 = M1 - 1;
    }
    // ClobVal <- ClobVal MIN TmpR;
    ClobVal = std::min(ClobVal, TmpR);

    // SAIL L408–L411: expand M2 right symmetrically
    // TmpR <- X[M2];
    // WHILE M2<I2 AND TmpR+Hyst GEQ X[M2+1] DO BEGIN TmpR <- TmpR MIN X[M2+1]; M2 <- M2+1 END;
    TmpR = X[M2];
    while (M2 < I2 && (TmpR + Hyst) >= X[M2 + 1]) {
      TmpR = std::min(TmpR, X[M2 + 1]);
      M2 = M2 + 1;
    }
    // ClobVal <- ClobVal MIN TmpR;
    ClobVal = std::min(ClobVal, TmpR);

    // SAIL L413: FOR i <- M1 STEP 1 UNTIL M2 DO X[i] <- ClobVal  (blank the peak region)
    for (std::int64_t i = M1; i <= M2; ++i) X[i] = ClobVal;

    // SAIL L416–L423: "Throw it back" — reject if peak is too narrow or at array boundary
    // IF (M2-M1+1 < MinWidth) OR MaxLoc=I1 OR MaxLoc=I2 THEN ...;
    if ((MinWidth != IntNoSet && (M2 - M1 + 1) < MinWidth) || (MaxLoc == I1) || (MaxLoc == I2)) {
      Nfound = Nfound - 1;
      continue;
    }
  }

  return static_cast<int>(Nfound);
}

} // namespace parshl
