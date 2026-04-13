// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "fftw_wrap.hpp"
#include "sail_array.hpp"
#include "parshl_dbffa.hpp"
#include "parshl_findpeaks.hpp"
#include <cmath>
#include <iostream>

int main() {
  const int N = 1024;

  parshl::FftwR2C fft(N);

  // Generate a sine
  const double fs = 48000.0;
  const double f0 = 1000.0;
  const double PI = 3.14159265358979323846;
  for (int n=0; n<N; ++n) {
    fft.in[n] = std::sin(2.0 * PI * f0 * (double)n / fs);
  }

  fft.exec();

  // Convert FFTW output to interleaved S[1:2*(N/2+1)]
  const long nBins = (N/2 + 1);
  sail::Array1D<double> S(1, 2*nBins);
  for (long i=1; i<=nBins; ++i) {
    S[2*i-1] = fft.out[i-1][0]; // Re
    S[2*i]   = fft.out[i-1][1]; // Im
  }

  sail::Array1D<double> Xdb(1, nBins);
  parshl::DbFromInterleavedComplex(S, 2*nBins, Xdb);

  sail::Array1D<double> Peaks(1, 20), Locs(1, 20);
  int nfound = parshl::FindPeaks(Xdb, Peaks, Locs, /*Thresh*/-300.0);
  std::cout << "Found " << nfound << " peaks\n";
  for (int i=1; i<=nfound; ++i) {
    std::cout << i << ": dB=" << Peaks[i] << " bin=" << Locs[i] << "\n";
  }
}
