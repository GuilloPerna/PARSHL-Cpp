// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "parshl_qisort.hpp"
#include <vector>
#include <algorithm>
//Archaeological note: SAIL uses quicksort, but the only semantic requirement is “ascending order” and carrying the second vector along. Stability does not harm fidelity; it only avoids arbitrary reorderings on ties.
namespace parshl {

void QIsort(sail::Array1D<double>& key, sail::Array1D<double>& val, std::int64_t n) {
  if (n <= 1) return;

  struct Pair { double k; double v; };
  std::vector<Pair> tmp;
  tmp.reserve(static_cast<std::size_t>(n));

  for (std::int64_t i = 1; i <= n; ++i) {
    tmp.push_back(Pair{ key[i], val[i] });
  }

  std::stable_sort(tmp.begin(), tmp.end(), [](const Pair& a, const Pair& b){
    return a.k < b.k;
  });

  for (std::int64_t i = 1; i <= n; ++i) {
    key[i] = tmp[static_cast<std::size_t>(i-1)].k;
    val[i] = tmp[static_cast<std::size_t>(i-1)].v;
  }
}

} // namespace parshl
