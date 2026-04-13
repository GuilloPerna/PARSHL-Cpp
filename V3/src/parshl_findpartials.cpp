// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port, archaeological extensions (V2), and experimental
//   improvements (V3): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "parshl_findpartials.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace parshl {

[[nodiscard]] static inline std::int64_t clamp_i64(std::int64_t v, std::int64_t lo, std::int64_t hi) noexcept {
  return std::max(lo, std::min(hi, v));
}

// SAIL source: INTEGER PROCEDURE FindPartials
// Parshl-source.txt L429–L497
//
// Wrapper around FindPeaks that adds frequency-range gating, minimum-separation
// enforcement, ascending-frequency sorting, and bin→Hz conversion.
// Called from the main analysis loop (Parshl-source.txt L1726–L1727).
//
// Steps mirroring SAIL:
//   L456: BinInt <- (MinSep/Fs)*Nfft + 0.5   (minimum separation in bins)
//   L458: ARRTRAN(Xsave, XmagDB)             (save/restore because FindPeaks blanks X)
//   L460–L461: I1 = 1+Nfft*(Fc1/Fs), I2 = 1+Nfft*(Fc2/Fs)  (band limits in bins)
//   L462: FindPeaks(...)                     (iterative maximum blanking)
//   L465: QIsort ascending by frequency      (needed so synthesiser can assign tracks)
//   L470: LinFrq[i] <- Fs*(LinFrq[i]-1)/Nfft  (bin→Hz conversion)
//   L463: ARRTRAN(XmagDB, Xsave)             (restore XmagDB for next call)
std::int64_t FindPartials(
  const sail::Array1D<double>& Xdb_in,
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
  bool trace
){
  // SAIL L458: ARRTRAN(Xsave, XmagDB) — save because FindPeaks modifies X in-place (blanking).
  // SAIL L463: ARRTRAN(XmagDB, Xsave) restores it; we use Xdb_in as const source instead.
  sail::Array1D<double> Xdb(Xdb_in.lo(), Xdb_in.hi());
  for (std::int64_t i = Xdb_in.lo(); i <= Xdb_in.hi(); ++i) Xdb[i] = Xdb_in[i];

  const std::int64_t Nspec = Xdb.hi();   // 1..Nspec  (SAIL L307: Nspec <- (Nfft DIV 2)+1)

  // SAIL L456: BinInt <- (MinSep/Fs)*Nfft + 0.5   (minimum peak separation in bins)
  const std::int64_t BinInt = static_cast<std::int64_t>((MinSep / Fs) * static_cast<double>(Nfft) + 0.5);

  // SAIL L460–L461: I1 = 1+Nfft*(Fc1/Fs),  I2 = 1+Nfft*(Fc2/Fs)  (search band, 1-based bins)
  std::int64_t I1 = 1 + static_cast<std::int64_t>(static_cast<double>(Nfft) * (Fc1 / Fs));
  std::int64_t I2 = 1 + static_cast<std::int64_t>(static_cast<double>(Nfft) * (Fc2 / Fs));
  I1 = clamp_i64(I1, 1, Nspec);
  I2 = clamp_i64(I2, 1, Nspec);
  if (I2 < I1) std::swap(I1, I2);

  // Work arrays (SAIL-style 1-based)
  sail::Array1D<double> Peaks(1, NpartialsReq);
  sail::Array1D<double> Locs(1, NpartialsReq);

  // SAIL L462: FindPeaks(Xdb, Peaks, Locs, ..., BinInt, I1, I2)
  // BlankWidth = BinInt enforces minimum separation between found peaks.
  const std::int64_t BlankWidth = BinInt;

  // SAIL L460: Npartials <- FindPeaks(...)  (REFERENCE INTEGER modified in-place in SAIL)
  std::int64_t Npartials = static_cast<std::int64_t>(parshl::FindPeaks(
    Xdb, Peaks, Locs,
    Thresh, Hyst,
    NpartialsReq,
    MinWid,
    BlankWidth,
    I1, I2
  ));

  // SAIL L465: QIsort(LinFrq, Npartials) — sort ascending by frequency (bin index).
  // Required so the synthesiser assigns oscillator tracks in spectral order.
  if (Npartials > 1) {
    parshl::QIsort(Locs, Peaks, Npartials);
  }

  // SAIL L470: LinFrq[Partial] <- Fs*(LinFrq[Partial]-1)/Nfft   (bin index → Hz)
  if (Npartials > NpartialsReq) Npartials = NpartialsReq;
  for (std::int64_t i = 1; i <= Npartials; ++i) {
    LinAmpDB[i] = Peaks[i];
    LinFrq[i] = Fs * ( (Locs[i] - 1.0) / (double)Nfft );
  }

  if (trace) {
    std::cout << "FindPartials: Npartials=" << Npartials
              << " I1=" << I1 << " I2=" << I2
              << " BinInt=" << BinInt << "\n";
  }

  return Npartials;
}

} // namespace parshl
