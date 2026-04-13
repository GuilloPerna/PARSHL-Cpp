// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#pragma once
#include <cstdint>
#include <vector>

namespace parshl {

struct TrackerParams {
  std::int64_t MaxOscs = 200;
  std::int64_t MaxLins = 200;

  // Runtime sample rate for forensic logging only (no decision logic should depend on this).
  double Fs = 0.0;

  // DFmax endpoints (PARSHL / SAIL)
  double Fc1 = 0.0;      // lower frequency limit
  double Fc2 = 22050.0;  // upper frequency limit (Nyquist typical)
  double DFmax1 = 0.0;   // max dev at low freq
  double DFmax2 = 0.0;   // max dev at high freq

  bool UDtrace = false;
  bool InstantRise = false;
  // V2: if true, disables oscillator recycling — exact SAIL monotone-Noscs behaviour.
  bool legacy_mode = false;
  // V2 D7: velocity-scale factor for adaptive DFmax (0.0 = disabled, D7 off).
  // DFmax_V2(osc) = max( DFmax_SAIL(osc),  dfmax_velocity_scale * |OscFrq[osc] - PrvOscFrq[osc]| )
  // PrvOscFrq is reset to OscFrq when PrvOscAmp==0 (OutAF rule), so new/recycled
  // oscillators automatically see velocity=0 and fall back to the SAIL static model.
  // Ignored (force 0) when legacy_mode = true.
  double dfmax_velocity_scale = 0.0;

  // V2 D6/DAmax: max tolerated amplitude change between frames, in dB (0.0 = disabled).
  // RECONSTRUCTIVE FEATURE — based on commented/disabled SAIL source:
  //   SAIL L158:  global declaration alongside DFmax: "REAL MinSep,...,DFmax1,DFmax2,DAmax,..."
  //   SAIL L962:  in Idefaults block: "# DAmax <- 10; # Disabled;"
  //   SAIL L1340 (JOS design note, in UpdateMap comment):
  //     "the previous amplitude tells whether this is a reasonable assumption [that the line
  //      disappeared], but the information is not being put to use here."
  // Semantic (§25/JOS-faithful): applied at no-line stop paths in get_closest_frq().
  //   Flag: PrvOscAmp[CurOsc] > 10^(-DAmax/20) at a no-line/no-match stop => damax_exceeded++.
  //   The oscillator ALWAYS stops (no line available); the counter is purely observational.
  //   This is the JOS insight: if an active oscillator loses its line, the previous amplitude
  //   tells whether the loss is acoustically surprising (loud) or expected (already quiet).
  // WARNING: units (dB) and threshold value are inferred, not confirmed by any SAIL execution
  // trace.  The SAIL value "10" was a disabled placeholder, not a calibrated default.
  // Recommended reference for experimental runs: 30 dB (10 dB counts only very loud stops).
  // NOTE: this comment block is superseded by the §26 DAmax_observe/DAmax_gate split below.

  // ---- D6 / DAmax: two independent variants (both OFF by default) ----
  //
  // Archaeological root:
  //   SAIL L158: global declaration: "REAL MinSep,...,DFmax1,DFmax2,DAmax,..."
  //   SAIL L962: commented, never executed: "# DAmax <- 10; # Disabled;"
  //   SAIL L1340 (JOS note in UpdateMap block comment):
  //     "the previous amplitude tells whether this is a reasonable assumption
  //      [that the line disappeared], but the information is not being put to use here."
  //
  // --- D6-observe (§25 / JOS-faithful) ---
  //   Applied at no-line stop paths ONLY (BestLin<=0, MinDist>DFmax).
  //   Does NOT alter the matching algorithm; does NOT cause extra stops.
  //   Purely observational: counts no-line stops where the oscillator was louder
  //   than the quiet floor 10^(-DAmax_observe/20).
  //   CLI: --damax-observe=N. Metric: damax_exceeded_observe.
  double DAmax_observe = 0.0;

  // --- D6-gate (§24 / reconstructive/experimental) ---
  //   Applied post-DFmax / pre-TakeIt (line was found and frequency-matched).
  //   Condition: |20*log10(LinAmp[BestLin] / PrvOscAmp[CurOsc])| > DAmax_gate => StopOsc.
  //   ALTERS the matching outcome: can stop oscillators that frequency-matching would accept.
  //   This is a reconstructive interpretation; NOT what the JOS comment describes.
  //   New/silent oscillators (PrvOscAmp==0) are always exempt.
  //   CLI: --damax-gate=N. Metric: damax_gate_stops.
  //   WARNING: not confirmed by any SAIL execution; units and insertion point are inferred.
  double DAmax_gate = 0.0;
  // [V3 E2] --closest-osc-free: when recycling dead oscillator slots, pick the slot
  // whose previous frequency (PrvOscFrq) is closest to the new line's frequency,
  // instead of the first available dead slot (V2 baseline).  Disabled by default;
  // V2 baseline is bit-identical when this is false.
  bool closest_osc_free = false;};

struct TrackerState {
  std::int64_t Noscs = 0;  // number of oscillators currently allocated

  // 1-based arrays sized Max+1
  std::vector<std::int64_t> LinOfOsc;   // >0 ON, 0 FREE, <0 SQUELCHED (stores -line)
  std::vector<std::int64_t> PrvLinOfOsc;
  std::vector<std::int64_t> OscOfLin;   // line -> osc (0 if free)

  std::vector<double> OscAmp, PrvOscAmp;
  std::vector<double> OscFrq, PrvOscFrq;

  // internal DFmax linear model (Alpha, Beta)
  bool df_ready = false;
  double alpha = 0.0;
  double beta = 0.0;
};

// Accumulated lifecycle / matching statistics.
// Non-algorithmic: purely observational counters reset in reset() and
// updated alongside the thread-local counters in update_map().
// Safe to read at any time; does not affect analysis or synthesis.
struct TrackerStats {
  std::int64_t add_count      = 0;  // free-to-ON or first-frame init transitions
  std::int64_t stop_count     = 0;  // ON-to-off (squelch/free) transitions
  std::int64_t recycle_count  = 0;  // V2: squelch-to-ON (recycled slot restart)
  std::int64_t dfmax_exceeded = 0;  // partials dropped: frequency jump > DFmax
  // D6-observe (§25/JOS-faithful): no-line stops where PrvOscAmp exceeded the quiet floor.
  // Purely observational — the stop is unconditional; no extra stops are introduced.
  // Always 0 when DAmax_observe == 0 (default/disabled).
  std::int64_t damax_exceeded_observe = 0;
  // D6-gate (§24/reconstructive): stops added by the post-DFmax amplitude gate.
  // These are ADDITIONAL stops beyond what frequency-matching alone would produce.
  // Always 0 when DAmax_gate == 0 (default/disabled).
  std::int64_t damax_gate_stops = 0;
  // [V3 E2] Number of recycled slots where closest-free slot differed from
  // first-free slot.  Always 0 when closest_osc_free == false (default).
  std::int64_t cof_reuse_count = 0;
};

class ParshlTracker {
public:
  explicit ParshlTracker(const TrackerParams& p);

  void reset();

  // Accumulated stats (non-algorithmic, observational).
  const TrackerStats& stats() const { return ts_; }

  // Frame1 is 1-based (if you have frame0, pass frame0+1)
  void update_map(
    std::int64_t Frame1,
    std::int64_t Nlins,
    const std::vector<double>& LinAmp, // lines amplitudes
    const std::vector<double>& LinFrq  // lines frequencies (Hz)
  );

  const TrackerState& state() const { return st_; }

  // [§32B] DBspec refinement API (JOS 8-JUL-85 conservative V2 approximation).
  // Used by --dbspec-refine in parshl_audio_peaks.cpp after update_map().
  // Read/write access to oscillator frequency only; no algorithmic side-effects.
  [[nodiscard]] bool osc_is_on(std::int64_t osc) const noexcept {
    const auto sz = static_cast<std::int64_t>(st_.LinOfOsc.size());
    if (osc < 1 || osc >= sz) return false;
    return st_.LinOfOsc[static_cast<std::size_t>(osc)] > 0;
  }
  [[nodiscard]] double osc_frq_hz(std::int64_t osc) const noexcept {
    const auto sz = static_cast<std::int64_t>(st_.OscFrq.size());
    if (osc < 1 || osc >= sz) return 0.0;
    return st_.OscFrq[static_cast<std::size_t>(osc)];
  }
  void set_osc_frq(std::int64_t osc, double frq_hz) noexcept {
    const auto sz = static_cast<std::int64_t>(st_.OscFrq.size());
    if (osc < 1 || osc >= sz) return;
    st_.OscFrq[static_cast<std::size_t>(osc)] = frq_hz;
  }
  // [V3 E3] Read/write access to oscillator amplitude; no algorithmic side-effects.
  // Used by --dbspec-refine-amp in parshl_audio_peaks.cpp after update_map().
  [[nodiscard]] double osc_amp_lin(std::int64_t osc) const noexcept {
    const auto sz = static_cast<std::int64_t>(st_.OscAmp.size());
    if (osc < 1 || osc >= sz) return 0.0;
    return st_.OscAmp[static_cast<std::size_t>(osc)];
  }
  void set_osc_amp(std::int64_t osc, double amp) noexcept {
    const auto sz = static_cast<std::int64_t>(st_.OscAmp.size());
    if (osc < 1 || osc >= sz) return;
    st_.OscAmp[static_cast<std::size_t>(osc)] = amp;
  }

  // [V4 T1] Inject an externally computed TrackerState.
  // Used by the peak-claiming tracker (parshl_peak_claim.cpp) to overwrite the
  // tracker state after running its own assignment logic, so that all downstream
  // synthesis and lifecycle-counting code can continue to read tracker.state()
  // without modification.
  // Contract: caller is responsible for a self-consistent TrackerState.
  //           tracker.stats() retains its previous value (V2 stats counters are
  //           not meaningful when T1 is active; T1Stats are reported separately).
  void overwrite_state(TrackerState st) noexcept { st_ = std::move(st); }

private:
  TrackerParams p_;
  TrackerState st_;
  TrackerStats ts_{};

  // ---- state predicates (match SAIL meaning) ----
  [[nodiscard]] bool prv_osc_on(std::int64_t osc) const noexcept { return st_.PrvLinOfOsc[osc] > 0; }
  [[nodiscard]] bool prv_osc_free(std::int64_t osc) const noexcept { return st_.PrvLinOfOsc[osc] == 0; }
  [[nodiscard]] bool prv_osc_squelched(std::int64_t osc) const noexcept { return st_.PrvLinOfOsc[osc] < 0; }

  [[nodiscard]] bool osc_on(std::int64_t osc) const noexcept { return st_.LinOfOsc[osc] > 0; }
  [[nodiscard]] bool osc_free(std::int64_t osc) const noexcept { return st_.LinOfOsc[osc] == 0; }
  [[nodiscard]] bool osc_squelched(std::int64_t osc) const noexcept { return st_.LinOfOsc[osc] < 0; }

  [[nodiscard]] bool lin_free(std::int64_t lin) const noexcept { return st_.OscOfLin[lin] == 0; }
  [[nodiscard]] bool lin_claimed(std::int64_t lin) const noexcept { return st_.OscOfLin[lin] > 0; }

  // ---- core helpers ----
  void stop_osc(std::int64_t osc);

  // DFmax according to SAIL CASE logic (Alpha/Beta computed once)
  double dfmax(std::int64_t osc);

  // robust accessor for LinFrq (handles 0-based vs 1-based vs "dummy slot")
  double lin_frq_at(std::int64_t lin, std::int64_t Nlins, const std::vector<double>& LinFrq) const;
  // NOTE: lin_amp_at is intentionally NOT declared here; D6/DAmax uses get_lin() directly
  // to avoid an unresolved-symbol trap (the method was declared but never defined in §23).
  // If a member-function form is ever needed, implement it fully in parshl_tracker.cpp first.

  // SAIL recursive matcher
  // D6-gate requires LinAmp (post-DFmax amplitude check).
  // D6-observe uses PrvOscAmp from tracker state; LinAmp is forwarded but only accessed
  // when DAmax_gate > 0. Passing it unconditionally keeps the signature stable.
  void get_closest_frq(
    std::int64_t CurOsc,
    double DFmax,
    std::int64_t Nlins,
    const std::vector<double>& LinFrq,
    const std::vector<double>& LinAmp
  );
  // V2: scan slots 1..Noscs for a dead (silent + not-ON) oscillator that can be recycled.
  // Returns slot index (1-based) or 0 if none found.
  std::int64_t find_free_osc() const noexcept;
  // [V3 E2] Like find_free_osc() but returns the dead slot whose PrvOscFrq is
  // closest to target_frq.  Slots with no frequency history (PrvOscFrq==0) are
  // treated as infinitely distant and lose to any slot with real history.
  // Falls back to find_free_osc() when all dead slots have no history.
  [[nodiscard]] std::int64_t find_closest_free_osc(double target_frq) const noexcept;
};

} // namespace parshl
