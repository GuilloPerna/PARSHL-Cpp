// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// parshl_t2_tracker.cpp
// [V4 T2] Global assignment tracker — prototype (§59); T2b trajectory cost (§60); T2c amplitude-rank cost (§61);
//         T2φ phase-aware cost (§66B).
//
// Origin:  Architecture A2 (docs/PARSHL_V4_ROADMAP.md §3 T2).
// Design:  docs/T2_GLOBAL_ASSIGNMENT_TRACKER.md (§58); docs/T2_PHASE_AWARE_COST.md (§66A).
// Impl:    §59; §60 (trajectory cost); §61 (amplitude-rank cost); §66B (phase-aware cost).
//
// See parshl_t2_tracker.hpp for the public API and contract.

#include "parshl_t2_tracker.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

namespace parshl {

// ─── Hungarian algorithm (O(N³)) ─────────────────────────────────────────────
//
// Solves the minimum-cost perfect matching on a square dense N×N cost matrix.
// Based on the classic dual-potential (Kuhn-Munkres) shortest-augmenting-path
// variant used in competitive programming (0-indexed interface).
//
// Input:
//   cost  — flat row-major N×N array (0-indexed); must be non-negative.
//   N     — matrix dimension (> 0).
//
// Returns:
//   ans[i] = column (0-indexed) assigned to row i, for i in [0..N).
//
// Notes:
//   - All values in `cost` must be finite (no ±∞ or NaN).
//   - LARGE/HUGE sentinels from t2_update() satisfy this; they are large finite
//     doubles, not IEEE infinities.
//   - Time complexity: O(N³).  Space: O(N²) for cost + O(N) auxiliary.
//
static std::vector<std::size_t> hungarian(
    const std::vector<double>& cost,
    std::size_t N)
{
  // Internally 1-based (rows/cols use [1..N]; slot [0] is a sentinel).
  const double INF = std::numeric_limits<double>::max() / 4.0;

  std::vector<double>      u(N + 1, 0.0);  // row potentials
  std::vector<double>      v(N + 1, 0.0);  // col potentials
  std::vector<std::size_t> p(N + 1, 0);    // p[j] = row assigned to col j
  std::vector<std::size_t> way(N + 1, 0);  // augmenting-path predecessor

  for (std::size_t i = 1; i <= N; ++i) {
    p[0] = i;               // sentinel: "trying to assign row i"
    std::size_t j0 = 0;
    std::vector<double> minv(N + 1, INF);
    std::vector<bool>   used(N + 1, false);

    // Shortest-augmenting-path inner loop.
    do {
      used[j0] = true;
      const std::size_t i0 = p[j0];
      std::size_t j1 = 0;
      double delta = INF;

      for (std::size_t j = 1; j <= N; ++j) {
        if (!used[j]) {
          // Reduced cost: c[i0][j] - u[i0] - v[j]
          const double cur = cost[(i0 - 1) * N + (j - 1)] - u[i0] - v[j];
          if (cur < minv[j]) {
            minv[j] = cur;
            way[j]  = j0;
          }
          if (minv[j] < delta) {
            delta = minv[j];
            j1    = j;
          }
        }
      }

      // Update potentials.
      for (std::size_t j = 0; j <= N; ++j) {
        if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
        else          { minv[j] -= delta; }
      }
      j0 = j1;
    } while (p[j0] != 0);

    // Trace back augmenting path and update assignment.
    do {
      const std::size_t j1 = way[j0];
      p[j0] = p[j1];
      j0    = j1;
    } while (j0);
  }

  // Convert to 0-based answer: ans[row] = col.
  std::vector<std::size_t> ans(N, 0);
  for (std::size_t j = 1; j <= N; ++j) {
    if (p[j] != 0) {
      ans[p[j] - 1] = j - 1;
    }
  }
  return ans;
}

// ─── t2_update ───────────────────────────────────────────────────────────────

T2Stats t2_update(
    const TrackerState&          prv,
    const TrackerParams&         tp,
    std::int64_t                 /*Frame1*/,
    std::int64_t                 Nlins,
    const std::vector<double>&   LinAmp,
    const std::vector<double>&   LinFrq,
    double                       wid,
    double                       wrank,
    std::vector<std::int64_t>&   last_peak_indices,
    double                       wphi,       // [§66B] phase cost weight (0 = disabled)
    double                       dt,         // [§66B] hop/fs
    const std::vector<double>&   LinPhase,   // [§66B] spectral phase per partial, 1-based
    std::vector<double>&         prv_phase,  // [§66B] per-slot previous spectral phase
    TrackerState&                out)
{
  T2Stats stats{};

  const auto MaxOscs    = tp.MaxOscs;
  const std::size_t sz  = static_cast<std::size_t>(MaxOscs + 1);   // 1-based
  const std::size_t lsz = static_cast<std::size_t>(tp.MaxLins + 1);

  // ── 1. Initialise out (same pattern as T1) ─────────────────────────────────
  out.Noscs    = 0;
  out.df_ready = prv.df_ready;
  out.alpha    = prv.alpha;
  out.beta     = prv.beta;

  out.LinOfOsc.assign(sz,  0LL);
  out.OscAmp  .assign(sz,  0.0);
  out.OscFrq  .assign(sz,  0.0);
  out.OscOfLin.assign(lsz, 0LL);  // not used by T2; kept zero

  out.PrvLinOfOsc.assign(sz, 0LL);
  out.PrvOscAmp  .assign(sz, 0.0);
  out.PrvOscFrq  .assign(sz, 0.0);

  const std::size_t prv_sz = std::min(sz, prv.LinOfOsc.size());
  for (std::size_t i = 0; i < prv_sz; ++i) {
    out.PrvLinOfOsc[i] = prv.LinOfOsc[i];
    if (i < prv.OscAmp.size()) out.PrvOscAmp[i] = prv.OscAmp[i];
    if (i < prv.OscFrq.size()) out.PrvOscFrq[i] = prv.OscFrq[i];
  }

  // [§60] Ensure last_peak_indices is sized for all slots (zero-extend for new slots).
  if (last_peak_indices.size() < sz) last_peak_indices.resize(sz, 0LL);
  // [§66B] Ensure prv_phase is sized for all slots (zero-extend for new slots).
  if (prv_phase.size() < sz) prv_phase.resize(sz, 0.0);

  // ── 2. Collect live oscillators ────────────────────────────────────────────
  std::vector<std::size_t> live_oscs;
  live_oscs.reserve(static_cast<std::size_t>(MaxOscs));
  for (std::size_t osc = 1; osc < prv_sz; ++osc) {
    if (prv.LinOfOsc[osc] > 0) live_oscs.push_back(osc);
  }

  const std::size_t M = live_oscs.size();
  const std::size_t K = (Nlins > 0) ? static_cast<std::size_t>(Nlins) : 0;

  // ── Helper: allocate a free oscillator slot for a new oscillator ───────────
  auto alloc_slot = [&](std::size_t pk) -> bool {
    if (pk < 1 || pk >= LinFrq.size()) return false;
    for (std::size_t osc = 1; osc < sz; ++osc) {
      const std::int64_t prv_lof = (osc < prv.LinOfOsc.size())
                                   ? prv.LinOfOsc[osc] : 0LL;
      if (out.LinOfOsc[osc] == 0 && prv_lof <= 0) {
        out.LinOfOsc[osc] = static_cast<std::int64_t>(pk);
        out.OscFrq[osc]   = LinFrq[pk];
        out.OscAmp[osc]   = LinAmp[pk];
        ++stats.creations;
        return true;
      }
    }
    return false;  // bank full; peak silently dropped
  };

  // ── Helper: recompute Noscs ────────────────────────────────────────────────
  auto compute_noscs = [&]() {
    std::int64_t noscs = 0;
    for (std::int64_t osc = MaxOscs; osc >= 1; --osc) {
      const std::size_t oi = static_cast<std::size_t>(osc);
      if (oi < out.LinOfOsc.size() && out.LinOfOsc[oi] > 0) {
        noscs = osc;
        break;
      }
    }
    out.Noscs = noscs;
  };

  // ── Edge cases ─────────────────────────────────────────────────────────────
  if (M == 0 && K == 0) { return stats; }

  if (K == 0) {
    // No peaks: terminate all live oscillators.
    stats.terminations += static_cast<std::int64_t>(M);
    // [§60] Reset trajectory memory for all live slots.
    for (std::size_t i = 0; i < M; ++i) {
      const std::size_t osc = live_oscs[i];
      if (osc < last_peak_indices.size()) last_peak_indices[osc] = 0LL;
    }
    return stats;
  }

  if (M == 0) {
    // No live oscillators: create one for each peak (up to bank limit).
    for (std::size_t pk = 1; pk <= K; ++pk) alloc_slot(pk);
    compute_noscs();
    return stats;
  }

  // ── [§61] Amplitude-rank penalty for each peak column ────────────────────────────
  //
  // rank_penalty[j] = rank(j) / max(1, K-1),  rank 0 = strongest peak.
  // Stronger peaks are cheaper to assign; weaker peaks carry extra cost.
  // With wrank=0 (flag off) all penalties stay 0.0 → §59/§60 behaviour unchanged.
  std::vector<double> rank_penalty(K, 0.0);
  if (wrank > 0.0 && K > 1) {
    std::vector<std::size_t> amp_order(K);
    for (std::size_t j = 0; j < K; ++j) amp_order[j] = j;
    std::sort(amp_order.begin(), amp_order.end(), [&](std::size_t a, std::size_t b) {
      const std::size_t pa = a + 1, pb = b + 1;
      const double aa = (pa < LinAmp.size()) ? LinAmp[pa] : 0.0;
      const double ab = (pb < LinAmp.size()) ? LinAmp[pb] : 0.0;
      return aa > ab;  // descending amplitude
    });
    const double denom = static_cast<double>(K - 1);
    for (std::size_t r = 0; r < K; ++r) {
      rank_penalty[amp_order[r]] = static_cast<double>(r) / denom;
    }
  }

  // ── 3. Build N×N cost matrix, N = max(M, K) ────────────────────────────────
  //
  // Layout:
  //   Rows 0..M-1  : live oscillators.
  //   Rows M..N-1  : dummy osc rows (only exist when K > M; padding = LARGE).
  //   Cols 0..K-1  : detected peaks.
  //   Cols K..N-1  : dummy peak cols (only exist when M > K; padding = LARGE).
  //
  // Cost convention:
  //   LARGE (1e9) : padding sentinel; also used as "no-match" for dummy rows/cols.
  //   HUGE  (2e9) : DFmax gate — pair is invalid; algorithm prefers LARGE (terminate)
  //                 over HUGE (bad match).
  //   Valid pair  : wf * |Δfrq_Hz| + wa * |Δamp_lin|  (small, << LARGE).
  //
  const std::size_t N     = (M >= K) ? M : K;
  const double      LARGE = 1.0e9;
  const double      HUGE  = 2.0e9;
  const double      wf    = 1.0;
  const double      wa    = 0.2;

  const double fc2 = (tp.Fc2 > 0.0) ? tp.Fc2 : 22050.0;
  auto dfmax_at = [&](double frq_hz) -> double {
    return tp.DFmax1 + (tp.DFmax2 - tp.DFmax1) * frq_hz / fc2;
  };

  // Allocate and fill with LARGE (padding / no-match default).
  std::vector<double> cost_mat(N * N, LARGE);

  for (std::size_t i = 0; i < M; ++i) {
    const std::size_t osc     = live_oscs[i];
    const double      osc_frq = (osc < prv.OscFrq.size()) ? prv.OscFrq[osc] : 0.0;
    const double      osc_amp = (osc < prv.OscAmp.size()) ? prv.OscAmp[osc] : 0.0;

    for (std::size_t j = 0; j < K; ++j) {
      const std::size_t pk = j + 1;  // 1-based peak index
      if (pk >= LinFrq.size()) continue;

      const double pk_frq = LinFrq[pk];
      const double pk_amp = LinAmp[pk];
      const double df     = std::abs(pk_frq - osc_frq);
      const double da     = std::abs(pk_amp - osc_amp);
      const double dfy    = dfmax_at(pk_frq);

      if (df <= dfy) {
        // [§60] Identity-switch penalty: add wid if this peak differs from last-known peak.
        const std::int64_t last_pk = (osc < last_peak_indices.size())
                                     ? last_peak_indices[osc] : 0LL;
        const double switch_pen = (wid > 0.0 && last_pk > 0
                                   && last_pk != static_cast<std::int64_t>(pk))
                                  ? wid : 0.0;
        // [§61] Amplitude-rank penalty: weaker peaks are slightly more expensive.
        const double rpen = (wrank > 0.0 && j < rank_penalty.size())
                            ? wrank * rank_penalty[j] : 0.0;
        // [§66B] Phase-coherence cost: |wrap(phi_obs - phi_predicted)|.
        // Gate (Strategy B): active only when oscillator has phase history
        //   (prv.PrvOscAmp[osc] > 0 => oscillator was live in previous frame).
        // phi_predicted = prv_phase[osc] + 2*pi * osc_frq * dt.
        // wrap result to [0, pi] via modulo 2pi.
        double phase_pen = 0.0;
        if (wphi > 0.0
            && osc < prv.PrvOscAmp.size() && prv.PrvOscAmp[osc] > 0.0
            && osc < prv_phase.size()
            && pk < LinPhase.size()) {
          constexpr double two_pi = 2.0 * std::numbers::pi;
          const double phi_obs  = LinPhase[pk];
          const double phi_pred = prv_phase[osc]
                                  + two_pi * osc_frq * dt;
          double phi_err = std::fmod(std::fabs(phi_obs - phi_pred), two_pi);
          if (phi_err > std::numbers::pi) phi_err = two_pi - phi_err;
          phase_pen = wphi * phi_err;
          ++stats.phase_eligible;
          stats.phase_cost_sum += phi_err;
        }
        cost_mat[i * N + j] = wf * df + wa * da + switch_pen + rpen + phase_pen;
      } else {
        cost_mat[i * N + j] = HUGE;
      }
    }
  }
  // Rows M..N-1 remain LARGE (dummy osc rows, all padding).
  // Cols K..N-1 remain LARGE (dummy peak cols, all padding).

  // ── 4. Solve LAP ──────────────────────────────────────────────────────────
  const auto asgn = hungarian(cost_mat, N);  // asgn[row] = col (0-based)

  // ── 5. Apply assignment ────────────────────────────────────────────────────
  std::vector<bool> peak_matched(K, false);

  for (std::size_t i = 0; i < M; ++i) {
    const std::size_t j   = asgn[i];        // assigned column (0-based)
    const std::size_t osc = live_oscs[i];

    if (j < K) {
      // Assigned to a real peak column.
      const double c = cost_mat[i * N + j];
      if (c < LARGE) {
        // Valid match (within DFmax).
        const std::size_t pk = j + 1;  // 1-based
        if (pk < LinFrq.size()) {
          out.LinOfOsc[osc] = static_cast<std::int64_t>(pk);
          out.OscFrq[osc]   = LinFrq[pk];
          out.OscAmp[osc]   = LinAmp[pk];
          peak_matched[j]   = true;
          ++stats.assignments;
          // [§60] Track identity switches and update last_peak_indices.
          const std::int64_t ipk = static_cast<std::int64_t>(pk);
          if (osc < last_peak_indices.size()) {
            if (last_peak_indices[osc] > 0 && last_peak_indices[osc] != ipk)
              ++stats.identity_switches;
            last_peak_indices[osc] = ipk;
          }
          // [§66B] Update stored spectral phase for next frame.
          if (osc < prv_phase.size() && pk < LinPhase.size())
            prv_phase[osc] = LinPhase[pk];
        } else {
          ++stats.terminations;  // bounds guard (should not occur)
        }
      } else {
        // DFmax-gated (cost == HUGE >= LARGE): osc terminates, peak freed.
        ++stats.terminations;
        if (osc < last_peak_indices.size()) last_peak_indices[osc] = 0LL;  // [§60]
        if (osc < prv_phase.size()) prv_phase[osc] = 0.0;                  // [§66B]
        // peak_matched[j] stays false → creation handled below.
      }
    } else {
      // Assigned to a padded (dummy peak) column: terminate oscillator.
      ++stats.terminations;
      if (osc < last_peak_indices.size()) last_peak_indices[osc] = 0LL;  // [§60]
      if (osc < prv_phase.size()) prv_phase[osc] = 0.0;                  // [§66B]
    }
  }

  // Create new oscillators for unmatched peaks.
  for (std::size_t j = 0; j < K; ++j) {
    if (!peak_matched[j]) alloc_slot(j + 1);
  }

  // ── 6. Recompute Noscs ────────────────────────────────────────────────────
  compute_noscs();
  return stats;
}

} // namespace parshl
