// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
#include "fftw_wrap.hpp"
#include "sail_array.hpp"
#include "parshl_getwin.hpp"
#include <vector>
#include <stdexcept>

namespace parshl {

// STFT frame runner aligned with PARSHL's separation of:
//   Nx = frame/window length (number of input samples per FFT frame)
//   Nfft = FFT length (power of two), with zero-padding if Nfft > Nx
//
// lastS is interleaved complex spectrum (1-based) for bins 1..(Nfft/2+1).
struct StftFrameRunner {
  int Nfft = 0;
  int Nx   = 0;
  int hop  = 0;

  std::vector<double> win; // length Nx

  FftwR2C fft;             // size Nfft
  sail::Array1D<double> lastS; // [1 .. 2*(Nfft/2+1)]

  explicit StftFrameRunner(
      int Nfft_,
      int hop_,
      int wintype = parshl::Hamming,
      int Nx_ = -1,
      double P3 = -1.0,
      double P4 = -1.0);

  [[nodiscard]] long nBins() const noexcept { return (Nfft/2 + 1); }

  // Compute dB magnitude spectrum into Xdb[1..nBins()], clobber-free.
  void compute_db(const std::vector<double>& x, long pos, sail::Array1D<double>& Xdb);
};

} // namespace parshl
