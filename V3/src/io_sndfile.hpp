// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port, archaeological extensions (V2), and experimental
//   improvements (V3): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
#include <string>
#include <vector>

namespace parshl {

// Read mono (or take channel 0) from any libsndfile-supported file.
// Output samples as double in [-1,1] (or whatever libsndfile provides), sample rate returned in fs.
bool ReadMonoAsDouble(const std::string& path, std::vector<double>& out, int& fs);

// Write mono float32 WAV.
bool WriteMonoWavFloat32(const std::string& path, const std::vector<double>& in, int fs);

} // namespace parshl
