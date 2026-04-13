// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port, archaeological extensions (V2), and experimental
//   improvements (V3): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// src/parshl_globals.hpp
#pragma once
#include <cstdint>

namespace parshl {

// ORIG: PARSHL SAIL defines:
//   IntNoSet = (1 LSH 35)  (sign-bit sentinel in 36-bit word)
//   RealNoSet = (-1.0@35)  (i.e., -1e35)
// We emulate the "negative sentinel" behavior explicitly in C++.
inline constexpr std::int64_t IntNoSet = -(std::int64_t(1) << 35);
inline constexpr double RealNoSet = -1.0e35;

// Debug flags (ORIG: Debug1/Debug3 globals used by FindPeaks/Qinterp)
extern bool Debug1;
extern bool Debug3;

} // namespace parshl
