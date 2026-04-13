// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port, archaeological extensions (V2), and experimental
//   improvements (V3): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "parshl_stft.hpp"
#include "parshl_dbffa.hpp"
#include <algorithm>

namespace parshl {

// SAIL source: per-analysis-loop initialisation before main "RI" loop
// Parshl-source.txt L1677: GetWin(WinBuf, WinType, Nx)
//
// Initialises the per-frame STFT runner:
//   - calls GetWin to build the analysis window of length Nx (L1677)
//   - pre-plans the FFTW r2c transform of size Nfft
//   - allocates lastS[1:2*(Nfft/2+1)] for interleaved Re/Im output
//     (matching the SAIL !FFA array layout, L1703)
StftFrameRunner::StftFrameRunner(int Nfft_, int hop_, int wintype, int Nx_, double P3, double P4)
: Nfft(Nfft_), Nx(Nx_ <= 0 ? Nfft_ : Nx_), hop(hop_),
  win(static_cast<size_t>(Nx_ <= 0 ? Nfft_ : Nx_), 1.0),
  fft(Nfft_),
  lastS(1, 2*(Nfft_/2 + 1))
{
  if (Nfft <= 0) throw std::runtime_error("StftFrameRunner: Nfft must be > 0");
  if (Nx <= 0)   throw std::runtime_error("StftFrameRunner: Nx must be > 0");
  if (Nx > Nfft) throw std::runtime_error("StftFrameRunner: Nx must be <= Nfft (zero-pad if needed)");

  // SAIL L1677: GetWin(WinBuf, WinType, Nx) — generate analysis window W[1:Nx]
  sail::Array1D<double> W(1, Nx);
  GetWin(W, wintype, Nx, P3, P4);
  for (int i = 1; i <= Nx; ++i) win[static_cast<size_t>(i-1)] = W[i];
}

// SAIL source: inner body of main analysis loop "RI"
// Parshl-source.txt L1695–L1717
//
// Per-frame STFT: zero-pad, window, FFT, compute XmagDB.
// Mirrors the SAIL sequence:
//   L1695: ARRCLR(X)                            — zero FFT input buffer
//   L1698: Sandi(InF(Channel), Xp, Nx, X, ...)  — read Nx samples into X[1:Nx]
//   L1701: IF WinType>1 THEN FOR i: X[i] <- X[i]*WinBuf[i]  — apply window
//   L1703: !FFA(X, Nfft)                         — in-place forward FFT (SIGLIB CFFT)
//   L1705–L1717: "DoFilter" loop                 — XmagDB = 10*log10(Re^2+Im^2  MAX  1e-20)
void StftFrameRunner::compute_db(const std::vector<double>& x, long pos, sail::Array1D<double>& Xdb) {
  // SAIL L1695: ARRCLR(X) — zero-initialise the Nfft-sample FFT input.
  // SAIL L1698/L1701: Sandi (read) then window multiplication X[i] <- X[i]*WinBuf[i].
  //
  // Zero-phase (ifftshift) layout — V1 CRITICAL correction:
  //   The window is centered at sample Mh=(Nx-1)/2.  Placing the windowed
  //   frame starting at buf[0] introduces a linear phase ramp of e^{-jπk(Nx-1)/Nfft}
  //   in every bin.  The SAIL/JOS zero-phase layout avoids this by splitting
  //   the frame around the FFT-buffer origin:
  //     positive half+centre: windowed[Mh..Nx-1] → buf[0..Mh]
  //     zeros in the middle:  buf[Mh+1..Nfft-Mh-1] = 0   (already set)
  //     negative half:        windowed[0..Mh-1]    → buf[Nfft-Mh..Nfft-1]
  for (int n = 0; n < Nfft; ++n) {
    fft.in[static_cast<size_t>(n)] = 0.0;
  }

  // Build windowed frame into a temporary buffer.
  std::vector<double> windowed(static_cast<size_t>(Nx), 0.0);
  for (int n = 0; n < Nx; ++n) {
    long idx = pos + n;
    double s = 0.0;
    if (idx >= 0 && idx < static_cast<long>(x.size())) s = x[static_cast<size_t>(idx)];
    windowed[static_cast<size_t>(n)] = s * win[static_cast<size_t>(n)];
  }

  const int Mh = (Nx - 1) / 2;

  // Positive half + centre: windowed[Mh..Nx-1] → buf[0..Mh]
  for (int n = 0; n <= Mh; ++n)
    fft.in[static_cast<size_t>(n)] = windowed[static_cast<size_t>(Mh + n)];

  // Negative half: windowed[0..Mh-1] → buf[Nfft-Mh..Nfft-1]
  for (int n = 0; n < Mh; ++n)
    fft.in[static_cast<size_t>(Nfft - Mh + n)] = windowed[static_cast<size_t>(n)];

  // SAIL L1703: !FFA(X, Nfft) — unnormalised forward FFT, in-place, interleaved Re/Im output.
  fft.exec();

  // SAIL L1703 output layout: X[2i-1] = Re(bin i),  X[2i] = Im(bin i),  i in {1..Nspec}
  // FFTW out[i-1][0..1] -> lastS[2i-1..2i]  (identical convention, 0-based vs 1-based).
  const long Nspec = (Nfft/2 + 1);   // SAIL L307: Nspec <- (Nfft DIV 2)+1
  for (long i = 1; i <= Nspec; ++i) {
    lastS[2*i - 1] = fft.out[static_cast<size_t>(i - 1)][0]; // Re
    lastS[2*i]     = fft.out[static_cast<size_t>(i - 1)][1]; // Im
  }

  // SAIL L1705–L1717: "DoFilter" — compute XmagDB from interleaved Re/Im.
  DbFromInterleavedComplex(lastS, 2*Nspec, Xdb);
}

} // namespace parshl
