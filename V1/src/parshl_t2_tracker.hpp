// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
// parshl_t2_tracker.hpp
// [V4 T2] Global assignment tracker — prototype (§59); T2b trajectory cost (§60); T2c amplitude-rank cost (§61);
//         T2φ phase-aware cost (§66B).
//
// Origin:  Architecture A2 (docs/PARSHL_V4_ROADMAP.md §3 T2).
// Design:  docs/T2_GLOBAL_ASSIGNMENT_TRACKER.md (§58); docs/T2_PHASE_AWARE_COST.md (§66A).
// Impl:    §59 (prototype); §60 (trajectory-aware cost); §61 (amplitude-rank cost); §66B (phase-aware cost).
//
// This module implements a global bipartite assignment tracker as an alternative
// to the PARSHL V2 greedy sequential tracker and the T1 peak-claiming tracker.
// Each frame is solved as a minimum-cost bipartite matching problem (LAP) via
// the classic O(N³) Hungarian algorithm.
//
// Contract:
//   - V2 baseline is bit-identical when --t2-global-assign is NOT present.
//   - Does NOT modify parshl_tracker.cpp, parshl_peak_claim.cpp, or any T1 code.
//   - Activated by --t2-global-assign in parshl_audio_peaks.cpp only.
//   - Not validated with --reverse-analysis or --t1-peak-claim.

#include "parshl_tracker.hpp"

#include <cstdint>
#include <vector>

namespace parshl {

// ---- Per-run accumulated statistics for the T2 global assignment tracker ----
// Emitted in [STATS_SUMMARY] as t2_assignments / t2_creations / t2_terminations.
struct T2Stats {
  std::int64_t assignments        = 0;  // existing oscillators matched to a peak (valid cost)
  std::int64_t creations          = 0;  // new oscillator slots allocated (peak without valid match)
  std::int64_t terminations       = 0;  // oscillators stopped (no peak matched within DFmax)
  std::int64_t identity_switches  = 0;  // [§60] assignments where peak identity changed from prv
  std::int64_t phase_eligible     = 0;  // [§66B] osc-peak pairs where phase term was active
  double       phase_cost_sum     = 0.0;// [§66B] sum of |phi_err| for all eligible pairs
};

// ---- t2_update — one frame of the global assignment tracker -----------------
//
// Replaces tracker.update_map() when --t2-global-assign is active.
//
// Algorithm (minimum-cost bipartite matching, O(N³)):
//   1. Collect live oscillators (count M) from prv.
//   2. Build N×N cost matrix, N = max(M, K):
//        cost(i, j) = wf * |Δfrq| + wa * |Δamp|                  if |Δfrq| <= DFmax(peak_j)
//                   + wid   * (j+1 != prv.LinOfOsc[osc])  [§60 trajectory term,  wid=0 → disabled]
//                   + wrank * rank_penalty(j)              [§61 amplitude-rank term, wrank=0 → disabled]
//        cost(i, j) = HUGE                                        otherwise (DFmax gate)
//        Padding cells (i >= M or j >= K) = LARGE
//   3. Solve LAP via Hungarian algorithm; obtain assignment[row] = col.
//   4. Apply assignment:
//        osc_i → real peak_j, cost < LARGE  →  valid match (assign).
//        osc_i → real peak_j, cost >= LARGE →  DFmax-gated: osc terminates, peak freed.
//        osc_i → padded col (j >= K)        →  osc terminates.
//        peak_j unmatched (peak_matched==false after step 4) → create new oscillator.
//   5. Recompute Noscs.
//
// Parameters:
//   prv    — TrackerState at the end of the previous frame (read-only).
//   tp     — TrackerParams (MaxOscs, MaxLins, DFmax1, DFmax2, Fc2 must be set).
//   Frame1 — 1-based frame index (same convention as update_map / t1_update).
//   Nlins  — number of partials in LinAmp / LinFrq (1-based: valid range [1..Nlins]).
//   LinAmp — partial amplitudes, 1-based [1..Nlins].
//   LinFrq — partial frequencies (Hz), 1-based [1..Nlins].
//   out    — output TrackerState; fully overwritten.
//
// Returns per-call T2Stats (caller accumulates across frames).
//
// Output state contract (same as T1):
//   out.PrvLinOfOsc[i] = prv.LinOfOsc[i]
//   out.PrvOscAmp[i]   = prv.OscAmp[i]
//   out.PrvOscFrq[i]   = prv.OscFrq[i]
//   out.LinOfOsc[i]    = peak-index if ON, 0 if FREE
//   out.OscFrq[i]      = assigned peak frequency if ON, 0 otherwise
//   out.OscAmp[i]      = assigned peak amplitude if ON, 0 otherwise
//   out.Noscs          = max active slot index
//
// Limitations (prototype):
//   - No oscillator recycling (no SQUELCH state); stopped slots return to FREE.
//   - dfmax_exceeded is not counted in tracker.stats() (V2 stats bypassed).
//   - No --closest-osc-free equivalent for new oscillator slot selection.
//   - Bank-full peaks (no free slot available) are silently dropped.
//   - Piano-C4 (N≈500): O(N³) per frame; expect ~30–90 s for full analysis.
// [§60] wid:   weight for identity-switch penalty in cost function (0.0 = disabled = §59 behaviour).
// [§61] wrank: weight for amplitude-rank penalty (0.0 = disabled = §60 behaviour).
//   rank_penalty(j) = rank(j) / max(1, K-1), rank 0 = highest-amplitude peak in current frame.
// last_peak_indices: persistent per-slot vector (sized MaxOscs+1) mapping slot → last assigned
//   peak index (1-based).  0 = new/unknown.  Updated by t2_update after each successful match.
//   Caller must preserve across frames; zero-initialised externally before first call.
// [§66B] wphi:  weight for phase-coherence cost term (0.0 = disabled = §61 behaviour).
//   c_phi(i,j) = |wrap(LinPhase[j] - (prv_phase[i] + 2π·osc_frq·dt))|, range [0,π].
//   Gate (Strategy B): phase term active only when prv.PrvOscAmp[osc] > 0 (oscillator has
//   phase history from previous matched frame).  New/recycled oscillators are exempt.
// [§66B] dt:   hop/fs (seconds per hop) — needed to compute expected phase advance.
// [§66B] LinPhase: spectral phase per partial (1-based), from atan2(Im, Re).  Size >= Nlins+1.
// [§66B] prv_phase: per-slot previous spectral phase (1-based, size MaxOscs+1).
//   Caller must preserve across frames.  Zero-initialised externally before first call.
//   Updated by t2_update: set to LinPhase[assigned_pk] on match; reset to 0.0 on terminate.
[[nodiscard]] T2Stats t2_update(
    const TrackerState&          prv,
    const TrackerParams&         tp,
    std::int64_t                 Frame1,
    std::int64_t                 Nlins,
    const std::vector<double>&   LinAmp,
    const std::vector<double>&   LinFrq,
    double                       wid,                // [§60] identity-switch penalty weight
    double                       wrank,              // [§61] amplitude-rank penalty weight
    std::vector<std::int64_t>&   last_peak_indices,  // [§60] per-slot last-assigned peak
    double                       wphi,               // [§66B] phase cost weight (0 = disabled)
    double                       dt,                 // [§66B] hop/fs (s)
    const std::vector<double>&   LinPhase,           // [§66B] spectral phase per partial, 1-based
    std::vector<double>&         prv_phase,          // [§66B] per-slot previous spectral phase
    TrackerState&                out
);

} // namespace parshl
