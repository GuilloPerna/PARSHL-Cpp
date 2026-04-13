// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// src/parshl_qinterp.hpp
#pragma once

namespace parshl {

// ORIG: INTERNAL SIMPLE REAL PROCEDURE Qinterp(REAL Ym1,Y0,Yp1; BOOLEAN InteriorX(TRUE))
// PARSHL listing: X <- (Yp1 - Ym1)/(2*(2*Y0 - Yp1 - Ym1)); clip to [-1,1] if InteriorX
[[nodiscard]] double Qinterp(double Ym1, double Y0, double Yp1, bool InteriorX = true);

} // namespace parshl
