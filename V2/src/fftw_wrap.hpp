// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
#include <fftw3.h>
#include <vector>
#include <stdexcept>

namespace parshl {

// Minimal wrapper for FFTW (double).
// No fancy RAII beyond ensuring plans are destroyed.
struct FftwR2C {
  int N = 0;
  fftw_plan plan = nullptr;
  std::vector<double> in;
  std::vector<fftw_complex> out;

  explicit FftwR2C(int n);
  ~FftwR2C();

  void exec();
};

struct FftwC2R {
  int N = 0;
  fftw_plan plan = nullptr;
  std::vector<fftw_complex> in;
  std::vector<double> out;

  explicit FftwC2R(int n);
  ~FftwC2R();

  void exec();
};

} // namespace parshl
