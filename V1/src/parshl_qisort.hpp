// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
#include <cstdint>
#include "sail_array.hpp"
//Archaeological note: SAIL uses quicksort, but the only semantic requirement is “ascending order” and carrying the second vector along. Stability does not harm fidelity; it only avoids arbitrary reorderings on ties.
namespace parshl {

// Sort key[1..n] ascending and carry val[1..n] alongside (stable).
void QIsort(sail::Array1D<double>& key, sail::Array1D<double>& val, std::int64_t n);

} // namespace parshl
