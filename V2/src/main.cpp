// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
//minimal compilable scaffold
#include "sail_array.hpp"
#include "parshl_findpeaks.hpp"
#include "parshl_globals.hpp"
#include <iostream>
#include "parshl_defaults.hpp"

auto Nfft  = parshl_defaults::Nfft;//this does not change logic; only ensures that if nothing is passed via CLI, the program starts the same as historical PARSHL.
auto Nx    = parshl_defaults::Nx;
auto Nhop  = parshl_defaults::Nhop;
auto Ndec  = parshl_defaults::Ndec;
auto Nh    = parshl_defaults::Nh;

auto winType = parshl_defaults::WinType;

bool doSynth     = parshl_defaults::DoSynth;
bool doFilter    = parshl_defaults::DoFlt;
bool swapOut     = parshl_defaults::SwapOut;
bool instantRise = parshl_defaults::InstantRise;

int traceLevel   = parshl_defaults::Trace;

int main() {
  std::cout << "PARSHL C++ port scaffold (archaeology mode)\n";

  // Mini smoke-test (no “improvements”; only to verify that it compiles and runs):
  // X[1:9] with a peak at the centre.
  sail::Array1D<double> X(1, 9);
  X[1]=0; X[2]=1; X[3]=2; X[4]=3; X[5]=10; X[6]=3; X[7]=2; X[8]=1; X[9]=0;

  sail::Array1D<double> Peaks(1, 10);
  sail::Array1D<double> Locs(1, 10);

  parshl::Debug1 = true;
  int n = parshl::FindPeaks(X, Peaks, Locs, /*Thresh*/0.0);

  std::cout << "Found " << n << " peak(s)\n";
  for (int i=1; i<=n; ++i) {
    std::cout << "  Peak " << i << ": val=" << Peaks[i] << " loc=" << Locs[i] << "\n";
  }
  return 0;
}
