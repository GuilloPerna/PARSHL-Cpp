// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once

#include <array>
#include <vector>

namespace parshl {

struct ParshlSynthState {
  static constexpr int SinSiz = 512;
  std::array<float, SinSiz> SinBuf{};
  double Mag = 0.0; // SinSiz/Fs
  bool sine_inited = false;
  std::vector<double> OscPhs; // 1-based: size MaxOscs+1
};

void parshl_synthesize_additive(
  ParshlSynthState& ss,
  int Nhop,
  int Bp,
  int Fs,
  int Frame1,
  int MaxOscs,
  const std::vector<double>& PrvOscAmp,
  const std::vector<double>& PrvOscFrq,
  const std::vector<double>& OscAmp,
  const std::vector<double>& OscFrq,
  int NskipActive,
  std::vector<float>& OutBuf,
  bool oscphs_double = false  // [V3 E4] if true: preserve fractional phase (no SAIL integer truncation)
);

} // namespace parshl
