// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port, archaeological extensions (V2), and experimental
//   improvements (V3): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "parshl_dbffa.hpp"
#include <cmath>
#include <algorithm>

namespace parshl {

// SAIL source: "DoFilter" block — XmagDB computation
// Parshl-source.txt L1686 (DBscl definition), L1705–L1717 ("DoFilter" loop body)
//
// Converts interleaved Re/Im FFT output (SAIL !FFA convention) to a dB power spectrum.
// Directly implements the XmagDB formula from the SAIL analysis loop:
//
//   L1686: DBscl <- (10.0/LOG(10.0))   [LOG = natural log in SIGLIB/SAIL]
//   L1710: XmagDB[k] <- DBscl * LOG(Xi*Xi + Xip1*Xip1  MAX  1.0@-20)
//
// Input S uses the SAIL/SIGLIB !FFA interleaved layout (Parshl-source.txt L1703):
//   S[2*k-1] = Re(bin k),  S[2*k] = Im(bin k),  k in {1..N/2}
//
// Output: Xdb[1..N/2] = 10*log10( max(Re(k)^2 + Im(k)^2,  1e-20) )
//
// FIX C (audited): noise floor = 1e-20  <-- from XmagDB analysis loop L1710
// DpyFFA (display-only, L310) used 1e-35 — that was the wrong source to copy.
void DbFromInterleavedComplex(
  const sail::Array1D<double>& S,
  long N,
  sail::Array1D<double>& Xdb
) {
  // SAIL L1686: DBscl <- (10.0/LOG(10.0))
  // SAIL note: "Incredibly, this is not performed at compile time"
  const double DBscl = 10.0 / std::log(10.0);

  // SAIL L1705–L1717: "DoFilter" loop
  //   XmagDB[(i+1)%2] <- DBscl * LOG(Xi*Xi + Xip1*Xip1  MAX  1.0@-20)
  // The (i+1)%2 indexing in SAIL maps interleaved pairs to output bins 1..Nspec.
  // Here we iterate directly by bin index (equivalent, just clearer C++ style).
  const long Nspec = N / 2;

  for (long I = 1; I <= Nspec; ++I) {
    const double Re = S[2*I - 1];
    const double Im = S[2*I];
    double pwr = Re*Re + Im*Im;
    if (pwr < 1.0e-20) pwr = 1.0e-20;
    Xdb[I] = DBscl * std::log(pwr);
  }
}

} // namespace parshl
