// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
// parshl_peak_claim.hpp
// [V4 T1] Peak-claiming tracker — alternative tracker for parshl-v4 research.
//
// Origin:  E5 failure analysis (§47C); Architecture A1 (docs/PARSHL_V3_FINAL.md §8).
// Design:  docs/T1_PEAK_CLAIMING_TRACKER.md (§54).
// Impl:    §55.
//
// This module implements a minimal peak-claiming tracker as an alternative to
// the PARSHL V2 greedy sequential tracker.  The assignment direction is inverted:
// peaks claim oscillators rather than oscillators searching for peaks.
//
// Contract:
//   - V2 baseline is bit-identical when --t1-peak-claim is NOT present.
//   - Does NOT modify parshl_tracker.cpp or any V2 code path.
//   - Activated by --t1-peak-claim in parshl_audio_peaks.cpp only.
//   - Not validated with --reverse-analysis (disabled when reverse_analysis=true).
//   - Not validated in combination with --dbspec-match, --grs-groups, or other V3 flags.

#include "parshl_tracker.hpp"

#include <cstdint>
#include <vector>

namespace parshl {

// ---- Per-run accumulated statistics for the T1 peak-claiming tracker ----------
// Emitted in [STATS_SUMMARY] as t1_created / t1_assigned / t1_terminated /
// t1_hold_survivals / t1_hold_terminations.
struct T1Stats {
  std::int64_t t1_created           = 0;  // new oscillator slots allocated (no prior match)
  std::int64_t t1_assigned          = 0;  // existing oscillators re-bound to a new peak
  std::int64_t t1_terminated        = 0;  // oscillators stopped (no peak claimed them)
  std::int64_t t1_hold_survivals    = 0;  // [§57] oscillators kept alive by squelch-hold
  std::int64_t t1_hold_terminations = 0;  // [§57] oscillators terminated after hold expired
};

// ---- t1_update — one frame of the peak-claiming tracker -----------------------
//
// Replaces tracker.update_map() when --t1-peak-claim is active.
//
// Algorithm (greedy, amplitude-priority):
//   1. Sort peaks by amplitude descending.
//   2. For each peak, find the closest ON oscillator (in prv) within DFmax Hz
//      that has not yet been claimed this frame.
//   3. Assign the oscillator to the peak (exclusive ownership).
//   4. If no candidate is available: allocate the first free slot.
//   5. ON oscillators in prv not claimed by any peak are terminated.
//
// Parameters:
//   prv    — TrackerState at the end of the previous frame (read-only).
//            On frame 1, this is the initial zero state from ParshlTracker::reset().
//   tp     — TrackerParams (MaxOscs, MaxLins, DFmax1, DFmax2, Fc2 must be set).
//   Frame1 — 1-based frame index (same convention as update_map).
//   Nlins  — number of partials in LinAmp / LinFrq (1-based: valid range [1..Nlins]).
//   LinAmp — partial amplitudes, 1-based [1..Nlins].
//   LinFrq — partial frequencies (Hz), 1-based [1..Nlins].
//   out    — output TrackerState; fully overwritten with the new state.
//
// Returns per-call T1Stats (caller accumulates across frames).
//
// Output state contract:
//   out.PrvLinOfOsc[i] = prv.LinOfOsc[i]    (lifecycle accounting preserved)
//   out.PrvOscAmp[i]   = prv.OscAmp[i]
//   out.PrvOscFrq[i]   = prv.OscFrq[i]
//   out.LinOfOsc[i]    = peak-index if ON, 0 if FREE  (no SQUELCH in prototype)
//   out.OscFrq[i]      = assigned peak frequency if ON, 0 otherwise
//   out.OscAmp[i]      = assigned peak amplitude if ON, 0 otherwise
//   out.Noscs          = max active slot index (= highest slot with LinOfOsc > 0)
//
// Limitations (prototype):
//   - No oscillator recycling (no SQUELCH state); stopped slots return to FREE.
//   - No --closest-osc-free equivalent.
//   - DFmax is the per-peak linear interpolation tp.DFmax1..tp.DFmax2 over [0..Fc2].
//   - Bank-full peaks (no free slot available) are silently dropped.
//
// [§57] squelch_hold and hold_counts:
//   squelch_hold = 0 → baseline T1 (immediate termination of unclaimed oscillators).
//   squelch_hold = N → unclaimed oscillators survive up to N consecutive frames before
//                      termination.  Sustained at their last freq/amp (constant hold).
//   hold_counts  — per-oscillator hold counter, 1-based, size MaxOscs+1.
//                  Must persist across frames; initialised/grown inside t1_update().
//                  Pass an empty vector on the first call; it will be resized as needed.
// [§64] t1_amp_hysteresis:
//   0.0  → no smoothing; A_new = A_peak (current behavior, bit-identical default).
//   0<N<1 → A_new = N*A_prev + (1-N)*A_peak  (exponential smoothing on matched slots).
//   Only applied when an oscillator is successfully matched to a peak.
//   Not applied on creation, termination, or squelch-hold frames.
[[nodiscard]] T1Stats t1_update(
    const TrackerState&          prv,
    const TrackerParams&         tp,
    std::int64_t                 Frame1,
    std::int64_t                 Nlins,
    const std::vector<double>&   LinAmp,
    const std::vector<double>&   LinFrq,
    int                          squelch_hold,  // [§57] 0=off, N=max hold frames
    std::vector<int>&            hold_counts,   // [§57] per-osc hold counter (persistent)
    double                       t1_amp_hysteresis,  // [§64] 0.0=off (default)
    TrackerState&                out
);

} // namespace parshl
