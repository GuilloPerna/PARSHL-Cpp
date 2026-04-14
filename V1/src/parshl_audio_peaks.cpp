// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// parshl_audio_peaks.cpp
// SAIL source: main analysis/synthesis driver — "RI" loop
// Parshl-source.txt L1690–L1780
//
// Replicates the SAIL main loop block-by-block:
//   L1695: ARRCLR(X) + read + window (via StftFrameRunner::compute_db)
//   L1703: !FFA(X, Nfft)             (via StftFrameRunner::compute_db)
//   L1705-L1717: XmagDB              (via DbFromInterleavedComplex)
//   L1726-L1727: FindPartials(...)   (via parshl::FindPartials)
//   L1728: LinAmp[i] <- 10^(dB/20)*SigScl  (inline below)
//   UpdateMap call                   (via tracker.update_map)
//   Synthesize                       (via parshl_synthesize_additive)
//   OLA flush IF Bp>=Nfft            (inline below)
//   L1945 (approx): SigScl <- 4/Nx  ("Guess for Hamming window, 50% overlap")
//
// Also contains:
//   decimate_sail(): SAIL INTEGER PROCEDURE Decimate (post-processing decimation)
//   InterpComplexAtPeak debug block: Paper §5 Step 4, NOT connected to pipeline
//
// Reads a WAV, runs a few frames and displays:
//
// - FindPartials (lines) + tracking
// - FindPeaks (peaks in dB) + locs
// - InterpComplexAtPeak (complex debug at bins)
//
// Args:
//   parshl_audio_peaks <input.wav> [Nfft=1024] [hop=256] [maxFrames=50] [thresh_db=-300]
//                    [MinSepHz=binWidth] [Hyst=0] [MinWid=1] [tracePartials=0] [Nx]
//   (Optional extra: to explicitly fix WinType in addition to Nx:
//      ... [tracePartials] [WinType] [Nx]   ==> only if passing *two* extra args)
//
// Examples:
//   ./parshl_audio_peaks flute-A5.wav 1024 256 1 -60
//   ./parshl_audio_peaks flute-A5.wav 1024 256 1 -60 43.066406 0 1 1 676

#include "io_sndfile.hpp"
#include "parshl_partials_writer.hpp"
#include "parshl_stft.hpp"
#include "parshl_findpeaks.hpp"
#include "parshl_complex_peak.hpp"
#include "parshl_findpartials.hpp"
#include "parshl_tracker.hpp"
#include "parshl_peak_claim.hpp"
#include "parshl_t2_tracker.hpp"
#include "parshl_synthesize.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <numbers>
#include <fstream>
#include <deque>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static void usage() {
  std::cerr
    << "usage:\n"
    << "  parshl_audio_peaks <input.wav> [Nfft=1024] [hop=256] [maxFrames=50] [thresh_db=-300]\n"
    << "                   [MinSepHz=binWidth] [Hyst=0] [MinWid=1] [tracePartials=0] [Nx]\n"
    << "  optional synth flags:\n"
    << "                   --legacy           (V1 mode: no complex interp, phase=0)\n"
    << "                   --synth [--flt --flt-ir <fir.txt>] [--synth-skip N] [--ndec N] --out-wav <path.wav> [--frames N]\n"
    << "  optional tracker DFmax (enables exploration without touching tracker algorithm):\n"
    << "                   --dfmax1 <Hz>  --dfmax2 <Hz>\n"
    << "  optional oscillator bank size (default 300):\n"
    << "                   --maxoscs=N        (N >= 10; isolates bank-capacity effect)\n"
    << "  optional explicit WinType:\n"
    << "                   ... [tracePartials] [WinType] [Nx]\n"
    << "  D4 reverse-time analysis (Serra & Smith 1990 §5, SAIL L1389-1421):\n"
    << "                   --reverse-analysis (analyse reversed audio, flip frame order for synthesis)\n"
    << "                   note: --reverse-analysis is incompatible with --legacy (ignored in legacy mode)\n"
    << "  D6 DAmax — two independent variants (both OFF by default, both suppressed by --legacy):\n"
    << "                   --damax-observe=N  (D6-observe: §25/JOS-faithful, purely observational, no extra stops)\n"
    << "                   --damax-gate=N     (D6-gate: §24/reconstructive, active post-DFmax gate, adds stops)\n"
    << "                   recommended: --damax-observe=30  --damax-gate=0 (gate is experimental)\n"
    << "  R1 OutPartials (Serra & Smith 1990 / SAIL Parshl-source.txt L500–510, L1818–1824):\n"
    << "                   --out-partials PREFIX  → writes PREFIX.amp.wav and PREFIX.frq.wav\n"
    << "                   float32 WAV, channels=MaxOscs, samplerate=round(Fs/hop) (frame rate)\n"
    << "                   channel k = oscillator k+1 amplitude (or frequency) at each frame\n"
    << "                   compatible with --reverse-analysis; does NOT affect tracker or synth\n"
    << "  R2 DBspec frequency refinement (JOS 8-JUL-85 / Parshl-source.txt L1425, conservative V2 approx):\n"
    << "                   --dbspec-refine    → post-match local max search in dB spectrum (frq only, no clobber)\n"
    << "                   note: matching unchanged; only OscFrq refined; no effect on amplitude, D4, or STATS when off\n"
    << "  E2 closest-osc-free (V3 experiment / Parshl-source.txt L1250-L1317 gap G10):\n"
    << "                   --closest-osc-free → when recycling dead slots, pick slot with PrvOscFrq closest to new line\n"
    << "                   (V2 baseline: first dead slot; this is experimental only)\n"
    << "  E1 SigScl analytic (V3 experiment / Parshl-source.txt ~L944 'IS SIGSCL CORRECT? ?'):\n"
    << "                   --sigscl-analytic  → SigScl = 2/sum(window) instead of SAIL heuristic 4/Nx\n"
    << "                   (DFT of real cosine: |X[k]| = (A/2)·Σw  →  SigScl=2/Σw; window-independent)\n"
    << "                   (only affects synthesis amplitude; tracker metrics unchanged)\n"
    << "  E4 OscPhs precision (V3 experiment / Parshl-source.txt L1046/L1667 audit C-2):\n"
    << "                   --oscphs-double    → preserve fractional phase at hop boundaries (no SAIL int truncation)\n"
    << "                   (only affects synthesis; tracker metrics unchanged)\n"
    << "  T1 peak-claiming tracker (V4 experiment / docs/T1_PEAK_CLAIMING_TRACKER.md §54):\n"
    << "                   --t1-peak-claim    → replace PARSHL greedy tracker with peak-claiming alternative\n"
    << "                   --t1-amp-hysteresis=N → [§64] exponential amplitude smoothing on matched oscillators\n"
    << "                                         N=0.0 (default) = off/bit-identical; N=0.2 = gentle; N=0.8 = heavy\n"
    << "                   (not validated with --reverse-analysis or other V3 flags)\n"
    << "  T2 global assignment tracker (V4 experiment / docs/T2_GLOBAL_ASSIGNMENT_TRACKER.md §58):\n"
    << "                   --t2-global-assign → replace greedy tracker with min-cost bipartite matching (Hungarian O(N³))\n"
    << "                   --t2-trajectory-cost → add wid*identity_switch penalty to Hungarian cost (auto-enables --t2-global-assign)\n"
    << "                   --t2-wid=N → set identity-switch penalty weight (default 5.0, requires --t2-trajectory-cost)\n"
    << "                   --t2-rank-cost → add amplitude-rank penalty to Hungarian cost (auto-enables --t2-global-assign)\n"
    << "                   --t2-wrank=N → set rank penalty weight (default 1.0, requires --t2-rank-cost)\n"
    << "                   (not validated with --reverse-analysis or --t1-peak-claim)\n"
    << "  N1 peak-normalize (post-synthesis, pre-write):\n"
    << "                   --normalize-peak=<dBFS>  → scale output so that |peak| == 10^(dBFS/20)\n"
    << "                   example: --normalize-peak=-1.0 prevents clipping on any synthesis\n"
    << "                   note: SDR is measured before normalization; does not affect analysis\n"
    << "\n"
    << "examples:\n"
    << "  parshl_audio_peaks flute-A5.wav 1024 256 12 -60\n"
    << "  parshl_audio_peaks flute-A5.wav 1024 256 1  -60 43.066406 0 1 1 676\n"
    << "  parshl_audio_peaks flute-A5.wav 1024 256 200 -60 --synth --reverse-analysis --out-wav out.wav\n";
}

[[nodiscard]] static inline std::int64_t BinInt_from_MinSepHz(double MinSepHz, double Fs, std::int64_t Nfft) {
  // PARSHL (SAIL) uses (MinSep/Fs)*Nfft + 0.5 with integer truncation.
  return (std::int64_t)((MinSepHz / Fs) * (double)Nfft + 0.5);
}

static inline bool trk_run_log_enabled() {
  const char* s = std::getenv("PARSHL_TRK_RUNLOG");
  return (s && *s && std::string(s) != "0");
}

static inline bool synth_log_enabled() {
  const char* s = std::getenv("PARSHL_SYNTH_LOG");
  return (s && *s && std::string(s) != "0");
}

// ---- Tracking + synthesis statistics accumulator --------------------------------
// Used to collect per-run metrics emitted in [STATS_SUMMARY] at the end of main.
// The struct is self-contained so it can be computed both from the live tracker
// (modes A/B) and from the pre-computed rev_frames sequence (mode C).
struct TrackAccum {
  // Partials / oscillator bank
  std::int64_t total_nlins    = 0;
  std::int64_t max_nlins      = 0;
  std::int64_t total_noscs    = 0;
  std::int64_t max_noscs      = 0;
  std::int64_t n_frames       = 0;
  // Lifecycle
  std::int64_t add_count      = 0;  // free→on transitions
  std::int64_t stop_count     = 0;  // on→off transitions
  std::int64_t recycle_count  = 0;  // squelch→on transitions
  std::vector<std::int64_t> lifetimes;  // per-trajectory duration in frames
  // Synthesis
  double peak_abs             = 0.0;
  std::int64_t overshoot_count = 0;
  double sum_sq               = 0.0;  // accumulates x^2 for RMS
  std::int64_t n_samples      = 0;
};

// Compute lifecycle metrics from rev_frames (mode C synthesis source).
// Detects transitions by comparing adjacent frames' LinOfOsc arrays.
// Still-alive trajectories at the end are closed with lifetime = end - birth + 1.
static TrackAccum compute_revframes_lifecycle(
    const std::vector<std::pair<std::int64_t, std::vector<std::int64_t>>>& frames_lof,
    // frames_lof[f] = {Noscs, LinOfOsc[0..MaxOscs]}
    std::int64_t MaxOscs,
    const std::vector<std::int64_t>& noscs_per_frame,
    const std::vector<std::int64_t>& nlins_per_frame)
{
  TrackAccum acc;
  acc.n_frames = static_cast<std::int64_t>(frames_lof.size());
  std::vector<std::int64_t> birth(static_cast<std::size_t>(MaxOscs + 1), -1);
  std::vector<std::int64_t> prev_lof(static_cast<std::size_t>(MaxOscs + 1), 0);

  for (std::size_t fi = 0; fi < frames_lof.size(); ++fi) {
    const auto& [noscs_cur, lof_cur] = frames_lof[fi];
    const std::int64_t nlins_cur = nlins_per_frame[fi];
    const std::int64_t f = static_cast<std::int64_t>(fi);

    acc.total_noscs += noscs_cur;
    acc.max_noscs = std::max(acc.max_noscs, noscs_cur);
    acc.total_nlins += nlins_cur;
    acc.max_nlins = std::max(acc.max_nlins, nlins_cur);

    const std::int64_t lim = std::min(MaxOscs, static_cast<std::int64_t>(lof_cur.size()) - 1);
    for (std::int64_t osc = 1; osc <= lim; ++osc) {
      const std::size_t oi = static_cast<std::size_t>(osc);
      const std::int64_t prv = prev_lof[oi];
      const std::int64_t cur = lof_cur[oi];
      if (prv <= 0 && cur > 0) {        // trajectory start
        if (prv < 0) ++acc.recycle_count;
        else         ++acc.add_count;
        birth[oi] = f;
      }
      if (prv > 0 && cur <= 0) {        // trajectory end
        ++acc.stop_count;
        if (birth[oi] >= 0) {
          acc.lifetimes.push_back(f - birth[oi]);
          birth[oi] = -1;
        }
      }
      prev_lof[oi] = cur;
    }
  }
  // Close trajectories still alive at the last frame.
  const std::int64_t last = acc.n_frames;
  for (std::int64_t osc = 1; osc <= MaxOscs; ++osc) {
    const std::size_t oi = static_cast<std::size_t>(osc);
    if (birth[oi] >= 0) {
      acc.lifetimes.push_back(last - birth[oi]);
    }
  }
  return acc;
}

static bool load_fir_coeffs(const std::string& path, std::vector<double>& ic) {
  ic.clear();
  std::ifstream in(path);
  if (!in) return false;
  std::string line;
  while (std::getline(in, line)) {
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    std::istringstream iss(line);
    double v = 0.0;
    while (iss >> v) ic.push_back(v);
  }
  return !ic.empty();
}

// SAIL source: INTEGER PROCEDURE Decimate
// Parshl-source.txt (utility procedure, post-processing path)
//
// Decimates array A[1:N] by factor M in-place, using stateful phase tracking
// (OWN INTEGER P preserves the fractional phase across successive calls).
// Returns the number of output samples written.
//
// SAIL listing:
//   INTEGER PROCEDURE Decimate(REAL ARRAY A; INTEGER N,M,I(0));
//     OWN INTEGER P;
//     IF M LEQ 1 THEN RETURN(N);
//     IF I NEQ 0 OR P LEQ 0 THEN P <- 1;
//     FOR i <- P STEP M UNTIL N DO A[j <- j+1] <- A[i];
//     P <- i-N;
//     RETURN(j);
[[nodiscard]] static int decimate_sail(std::vector<float>& A, int N, int M, int I, int& P) {
  if (M <= 1) return N;
  if (I != 0 || P <= 0) P = 1;

  int j = 0;
  int i = P;
  for (; i <= N; i += M) {
    A[static_cast<std::size_t>(j)] = A[static_cast<std::size_t>(i - 1)];
    ++j;
  }
  P = i - N;
  return j;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    usage();
    return 1;
  }

  const std::string inPath = argv[1];

  bool synth_enabled = false;
  bool flt_enabled = false;
  int synth_skip = 0;
  int ndec = 1;
  std::string flt_ir_path;
  std::string synth_out_wav;
  std::string out_partials_prefix; // --out-partials PREFIX → PREFIX.amp.wav + PREFIX.frq.wav
  std::int64_t frames_override = -1;
  // DFmax1/DFmax2 CLI overrides; negative sentinel = use driver default.
  double cli_dfmax1 = -1.0;
  double cli_dfmax2 = -1.0;
  // V2 D7: velocity-scale for adaptive DFmax (0.0 = disabled).
  double cli_dfmax_vscale = 0.0;
  // V2 mode by default; --legacy restores exact V1 (SAIL) behaviour:
  // integer-bin amplitude, phase=0 for all new oscillators.
  bool legacy_mode = false;
  // V2: oscillator bank capacity.  -1 = use default (300).
  // --maxoscs=N isolates the effect of bank size without touching the matching algorithm.
  std::int64_t cli_maxoscs = -1;
  // D6-observe (§25/JOS-faithful): counts no-line stops on loud oscillators. No extra stops.
  // CLI: --damax-observe=N  (default 0 = disabled). Metric: damax_exceeded_observe.
  double cli_damax_observe = 0.0;
  // D6-gate (§24/reconstructive): active amplitude gate post-DFmax/pre-TakeIt. Adds stops.
  // CLI: --damax-gate=N  (default 0 = disabled). Metric: damax_gate_stops.
  // WARNING: experimental; not confirmed by any SAIL execution trace.
  double cli_damax_gate = 0.0;
  // D4 reverse-time analysis (Serra & Smith 1990 §5 / SAIL L1389-1421).
  // When true: audio is reversed, full analysis pass is run on the reversed signal,
  // the resulting per-frame parameters are collected then flipped back to forward order.
  // Synthesis (if --synth) uses those forward-ordered parameters, giving the tracker
  // a "look-ahead" through the reversed audio so that attack transients are analysed
  // from their stable sustain end first.
  // Disabled when --legacy is active (--legacy takes priority).
  bool reverse_analysis = false;
  // [§32B] DBspec frequency refinement (JOS 8-JUL-85 conservative V2 approximation).
  // After update_map(), each active osc searches Xdb locally for a better frq estimate.
  // Xdb is read-only; matching/amplitude/D4/STATS are unaffected when disabled.
  bool dbspec_refine = false;
  // [V3 E2] --closest-osc-free: when recycling dead oscillator slots, pick the slot
  // whose PrvOscFrq is closest to the new line, instead of the first dead slot (V2 baseline).
  bool closest_osc_free = false;
  // [V3 E1] --sigscl-analytic: replace the SAIL heuristic SigScl=4/Nx with an
  // analytically derived value SigScl=2/sum(w[n]) based on the actual window.
  // Origin Type C: derived from DFT reconstruction theory.
  // Source A: Parshl-source.txt ~L944: "IS SIGSCL CORRECT? ? (Synth scaling)"
  // Source B: DFT of real cosine A·cos(2πk₀n/N) with window w:
  //           |X[k₀]| = (A/2)·Σw  →  LinAmp = |X[k₀]|·SigScl = A  →  SigScl = 2/Σw.
  //           Factor 2 is a DFT property (real signal); Σw adapts to the actual window.
  bool sigscl_analytic = false;
  // [V3 E3] --dbspec-refine-amp: post-match amplitude refinement from same Xdb local peak.
  // Origin Type B+C:
  //   B: JOS 8-JUL-85 (Parshl-source.txt L1425): "Consider allowing Oscs to find their peaks
  //      in DBspec instead of LinFrqs." — conservative V2 approximation.
  //   C: Analytical deduction: the Qinterp-interpolated dB peak directly yields amplitude
  //      via 10^(db/20)*SigScl (same conversion used in the LinAmp pipeline).
  // Contract: requires/auto-promotes --dbspec-refine; flag off = V2 baseline unchanged.
  bool dbspec_refine_amp = false;
  // [V3 E4] --oscphs-double: remove SAIL per-hop integer truncation of OscPhs.
  // Origin Type C: numerical precision — per-hop truncation discards fractional
  // phase (up to 1 table step = 1/512 cycle ≈ 0.7°) at each hop boundary.
  // Source: Parshl-source.txt L1046/L1667 INTEGER ARRAY (audit C-2).
  // Contract: experimental; flag off = SAIL-faithful truncation; tracker unchanged.
  bool oscphs_double = false;
  // [V3 E5] --dbspec-match: post-update_map candidate override via DBspec peak search.
  // Origin Type B+C:
  //   B: JOS 8-JUL-85 (Parshl-source.txt L1425–L1432):
  //      "Consider allowing Oscs to find their peaks in DBspec instead of LinFrqs."
  //   C: Conservative Variant A — maintain FindPartials+LinFrq+update_map;
  //      search Xdb within ±DFmax_bins around PrvOscFrq (predicted frequency);
  //      replace the LinFrq match only if the DBspec peak is strictly closer.
  // Metric: dbspec_match_count = total oscillator-frame replacements.
  // Contract: Experimental. Disabled by default. V2 baseline unchanged when off.
  bool dbspec_match = false;
  std::int64_t dbspec_match_count = 0;
  // [V3 E6] --grs-groups=N: GRS spectral group reduction.
  // Origin Type B: GRS comment (Parshl-source.txt L1432–L1436):
  //   "GRS suggested doing groups reduction on the mag² spectrum,
  //    letting each group behave as a 'virtual line'."
  // Algorithm: divide [Fc1,Fc2] into N equal-bandwidth groups; per group:
  //   centroid_frq = sum(frq_b * mag2_b) / sum(mag2_b);
  //   virtual_amp   = sqrt(sum_mag2) * SigScl.
  //   Replaces LinFrq_v[1..Nlins_eff]/LinAmp[1..Nlins_eff] for tracker.update_map.
  //   FindPartials still runs: Nlins reported in stats; LinPhase zeroed for virtual lines.
  // Metric: grs_group_collision_count — cumulative count of (frame×group) pairs
  //   where ≥2 FindPartials peaks fall in the same group (spectral ambiguity indicator).
  // Contract: Experimental. Disabled by default. V2 baseline unchanged when off.
  std::int64_t grs_groups = 0;
  std::int64_t grs_group_collision_count = 0;
  // [V4 T1] --t1-peak-claim: replace PARSHL greedy tracker with peak-claiming alternative.
  // Origin: E5 failure analysis (\u00a747C); Architecture A1 (docs/PARSHL_V3_FINAL.md \u00a78).
  // Design: docs/T1_PEAK_CLAIMING_TRACKER.md (\u00a754). Implementation: \u00a755.
  // Contract: Experimental. Disabled by default. V2 baseline unchanged when off.
  //           Not validated with --reverse-analysis or other V3 flags.
  bool t1_peak_claim = false;
  parshl::T1Stats t1_stats_accum{};
  int t1_squelch_hold = 0;       // [§57] 0 = off, N = max hold frames
  std::vector<int> t1_hold_counts;  // [§57] per-osc hold counter, persisted across frames
  double t1_amp_hysteresis = 0.0;   // [§64] 0=off (default), 0<N<1 = smoothing
  // [§65] Analysis error diagnostics — pure instrumentation, no synthesis change.
  // freq_error:  |f_interpolated - f_integer_bin| per partial (Hz) — Qinterp correction magnitude.
  // amp_error:   |LinAmp[p] - OscAmp[osc]| per matched oscillator (linear) — post-tracking drift.
  // phase_error: |phi_t - (phi_{t-1} + omega*dt)| per matched oscillator (rad, wrapped to [0,pi]).
  //              McAulay-Quatieri spectral phase continuity; needs no synthesis.
  double diag_freq_err_sum = 0.0;
  double diag_freq_err_max = 0.0;
  std::int64_t diag_freq_err_n = 0;
  double diag_amp_err_sum  = 0.0;
  std::int64_t diag_amp_err_n  = 0;
  double diag_phase_err_sum = 0.0;
  double diag_phase_err_max = 0.0;
  std::int64_t diag_phase_err_n = 0;
  std::vector<double> diag_prv_phase; // stored spectral phase per oscillator (1-based), previous frame
  // [V4 T2] --t2-global-assign: global min-cost assignment tracker.
  // Design: docs/T2_GLOBAL_ASSIGNMENT_TRACKER.md (§58). Implementation: §59.
  // Contract: Experimental. Disabled by default. V2 baseline unchanged when off.
  //           Not validated with --reverse-analysis or --t1-peak-claim.
  bool t2_global_assign = false;
  parshl::T2Stats t2_stats_accum{};
  bool t2_trajectory_cost = false;               // [§60] identity-switch penalty enabled
  double t2_wid = 5.0;                           // [§60] penalty weight (tunable via --t2-wid=)
  std::vector<std::int64_t> t2_last_peak_indices; // [§60] per-slot trajectory state
  bool   t2_rank_cost = false;                    // [§61] amplitude-rank penalty enabled
  double t2_wrank     = 1.0;                      // [§61] rank penalty weight (tunable via --t2-wrank=)
  bool   t2_phase_cost = false;                   // [§66B] phase-coherence cost term enabled
  double t2_wphi      = 1.0;                      // [§66B] phase cost weight (tunable via --t2-wphi=)
  std::vector<double> t2_prv_phase;               // [§66B] per-slot previous spectral phase
  // [§68] T4 phase-coherent synthesis
  bool   t4_phase_coherent = false;          // [§68] T4 enabled (S3 + S2)
  double t4_alpha          = 0.0;            // [§68] S2 correction weight in (0,1]
  std::int64_t t4_correction_count_accum = 0; // [§68] accumulated correction events
  // [N1] --normalize-peak=<dBFS>: post-synthesis peak-level normalization.
  // Scales the output buffer so that |peak| = 10^(dBFS/20) before writing the WAV.
  // Applied after the 32768 rescale; does NOT affect SDR measurement or tracker.
  bool   normalize_peak_enabled = false;
  double normalize_peak_dbfs    = -1.0;
  std::vector<std::string> pos_args;
  pos_args.reserve(static_cast<std::size_t>(std::max(0, argc - 2)));

  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--synth") {
      synth_enabled = true;
      continue;
    }
    if (a == "--flt") {
      flt_enabled = true;
      continue;
    }
    if (a == "--flt-ir") {
      if (i + 1 >= argc) {
        std::cerr << "missing value for --flt-ir\n";
        usage();
        return 1;
      }
      flt_ir_path = argv[++i];
      continue;
    }
    if (a.rfind("--flt-ir=", 0) == 0) {
      flt_ir_path = a.substr(std::string("--flt-ir=").size());
      continue;
    }
    if (a == "--synth-skip") {
      if (i + 1 >= argc) {
        std::cerr << "missing value for --synth-skip\n";
        usage();
        return 1;
      }
      synth_skip = std::stoi(argv[++i]);
      continue;
    }
    if (a.rfind("--synth-skip=", 0) == 0) {
      synth_skip = std::stoi(a.substr(std::string("--synth-skip=").size()));
      continue;
    }
    if (a == "--ndec") {
      if (i + 1 >= argc) {
        std::cerr << "missing value for --ndec\n";
        usage();
        return 1;
      }
      ndec = std::stoi(argv[++i]);
      continue;
    }
    if (a.rfind("--ndec=", 0) == 0) {
      ndec = std::stoi(a.substr(std::string("--ndec=").size()));
      continue;
    }
    if (a == "--out-wav") {
      if (i + 1 >= argc) {
        std::cerr << "missing value for --out-wav\n";
        usage();
        return 1;
      }
      synth_out_wav = argv[++i];
      continue;
    }
    if (a.rfind("--out-wav=", 0) == 0) {
      synth_out_wav = a.substr(std::string("--out-wav=").size());
      continue;
    }
    if (a == "--frames") {
      if (i + 1 >= argc) {
        std::cerr << "missing value for --frames\n";
        usage();
        return 1;
      }
      frames_override = std::stoll(argv[++i]);
      continue;
    }
    if (a.rfind("--frames=", 0) == 0) {
      frames_override = std::stoll(a.substr(std::string("--frames=").size()));
      continue;
    }
    if (a == "--dfmax1") {
      if (i + 1 >= argc) { std::cerr << "missing value for --dfmax1\n"; usage(); return 1; }
      cli_dfmax1 = std::stod(argv[++i]);
      continue;
    }
    if (a.rfind("--dfmax1=", 0) == 0) {
      cli_dfmax1 = std::stod(a.substr(std::string("--dfmax1=").size()));
      continue;
    }
    if (a == "--dfmax2") {
      if (i + 1 >= argc) { std::cerr << "missing value for --dfmax2\n"; usage(); return 1; }
      cli_dfmax2 = std::stod(argv[++i]);
      continue;
    }
    if (a.rfind("--dfmax2=", 0) == 0) {
      cli_dfmax2 = std::stod(a.substr(std::string("--dfmax2=").size()));
      continue;
    }
    if (a == "--legacy") {
      legacy_mode = true;
      continue;
    }
    if (a == "--reverse-analysis") {
      reverse_analysis = true;
      continue;
    }
    if (a == "--dfmax-vscale") {
      if (i + 1 >= argc) { std::cerr << "missing value for --dfmax-vscale\n"; usage(); return 1; }
      cli_dfmax_vscale = std::stod(argv[++i]);
      continue;
    }
    if (a.rfind("--dfmax-vscale=", 0) == 0) {
      cli_dfmax_vscale = std::stod(a.substr(std::string("--dfmax-vscale=").size()));
      continue;
    }
    // D6-observe (§25/JOS-faithful): quiet-floor counter at no-line stops. No extra stops.
    if (a == "--damax-observe") {
      if (i + 1 >= argc) { std::cerr << "missing value for --damax-observe\n"; usage(); return 1; }
      cli_damax_observe = std::stod(argv[++i]);
      continue;
    }
    if (a.rfind("--damax-observe=", 0) == 0) {
      cli_damax_observe = std::stod(a.substr(std::string("--damax-observe=").size()));
      continue;
    }
    // D6-gate (§24/reconstructive): active amplitude gate post-DFmax/pre-TakeIt. Adds stops.
    if (a == "--damax-gate") {
      if (i + 1 >= argc) { std::cerr << "missing value for --damax-gate\n"; usage(); return 1; }
      cli_damax_gate = std::stod(argv[++i]);
      continue;
    }
    if (a.rfind("--damax-gate=", 0) == 0) {
      cli_damax_gate = std::stod(a.substr(std::string("--damax-gate=").size()));
      continue;
    }
    // V2: configurable oscillator bank size — isolates bank-capacity effect.
    if (a == "--maxoscs") {
      if (i + 1 >= argc) { std::cerr << "missing value for --maxoscs\n"; usage(); return 1; }
      cli_maxoscs = std::stoll(argv[++i]);
      continue;
    }
    if (a.rfind("--maxoscs=", 0) == 0) {
      cli_maxoscs = std::stoll(a.substr(std::string("--maxoscs=").size()));
      continue;
    }
    // R1 OutPartials — §28: write per-frame OscAmp/OscFrq to PREFIX.amp.wav / PREFIX.frq.wav.
    // No effect on tracker algorithm, synthesis, or any console output.
    if (a == "--out-partials") {
      if (i + 1 >= argc) {
        std::cerr << "missing value for --out-partials\n";
        usage();
        return 1;
      }
      out_partials_prefix = argv[++i];
      continue;
    }
    if (a.rfind("--out-partials=", 0) == 0) {
      out_partials_prefix = a.substr(std::string("--out-partials=").size());
      continue;
    }
    // [§32B] DBspec frequency refinement (JOS 8-JUL-85 conservative V2 approximation).
    // Parshl-source.txt L1425: "Consider allowing Oscs to find their peaks in DBspec
    // instead of LinFrqs." -- NOT the full DBspec matching; only post-match freq refinement.
    if (a == "--dbspec-refine") {
      dbspec_refine = true;
      continue;
    }
    // [V3 E2] Closest-free-oscillator dispatch.
    if (a == "--closest-osc-free") {
      closest_osc_free = true;
      continue;
    }
    // [V3 E1] SigScl analytic normalisation.
    if (a == "--sigscl-analytic") {
      sigscl_analytic = true;
      continue;
    }
    // [V3 E3] DBspec amplitude refinement (post-match OscAmp update from Xdb peak).
    if (a == "--dbspec-refine-amp") {
      dbspec_refine_amp = true;
      continue;
    }
    // [V3 E4] OscPhs precision — remove SAIL integer truncation at hop boundaries.
    if (a == "--oscphs-double") {
      oscphs_double = true;
      continue;
    }
    // [V3 E5] DBspec matching: post-update_map candidate override via Xdb peak search.
    if (a == "--dbspec-match") {
      dbspec_match = true;
      continue;
    }
    // [V3 E6] GRS spectral group reduction.
    if (a.rfind("--grs-groups=", 0) == 0) {
      grs_groups = std::stoll(a.substr(std::string("--grs-groups=").size()));
      continue;
    }
    if (a == "--grs-groups") {
      if (i + 1 >= argc) { std::cerr << "missing value for --grs-groups\n"; usage(); return 1; }
      grs_groups = std::stoll(argv[++i]);
      continue;
    }
    // [V4 T1] Peak-claiming tracker.
    if (a == "--t1-peak-claim") {
      t1_peak_claim = true;
      continue;
    }
    if (a.rfind("--t1-amp-hysteresis=", 0) == 0) {
      t1_amp_hysteresis = std::stod(a.substr(std::string("--t1-amp-hysteresis=").size()));
      t1_peak_claim = true;  // auto-enable T1 when hysteresis is set
      continue;
    }
    if (a == "--t1-amp-hysteresis") {
      if (i + 1 >= argc) { std::cerr << "missing value for --t1-amp-hysteresis\n"; return 1; }
      t1_amp_hysteresis = std::stod(argv[++i]);
      t1_peak_claim = true;
      continue;
    }
    // [V4 T1 §57] Squelch-hold continuity patch.
    if (a.rfind("--t1-squelch-hold=", 0) == 0) {
      t1_squelch_hold = std::stoi(a.substr(std::string("--t1-squelch-hold=").size()));
      continue;
    }
    if (a == "--t1-squelch-hold") {
      if (i + 1 >= argc) { std::cerr << "missing value for --t1-squelch-hold\n"; usage(); return 1; }
      t1_squelch_hold = std::stoi(argv[++i]);
      continue;
    }
    // [V4 T2] Global assignment tracker.
    if (a == "--t2-global-assign") {
      t2_global_assign = true;
      continue;
    }
    // [§60] Trajectory-cost extension: identity-switch penalty.
    if (a == "--t2-trajectory-cost") {
      t2_trajectory_cost = true;
      t2_global_assign   = true;  // auto-enable the parent feature
      continue;
    }
    if (a.rfind("--t2-wid=", 0) == 0) {
      t2_wid = std::stod(a.substr(9));
      continue;
    }
    // [§61] Amplitude-rank cost extension.
    if (a == "--t2-rank-cost") { t2_rank_cost = true; t2_global_assign = true; continue; }
    if (a.rfind("--t2-wrank=", 0) == 0) { t2_wrank = std::stod(a.substr(11)); continue; }
    // [§66B] Phase-coherence cost extension.
    if (a == "--t2-phase-cost") { t2_phase_cost = true; t2_global_assign = true; continue; }
    if (a.rfind("--t2-wphi=", 0) == 0) { t2_wphi = std::stod(a.substr(10)); continue; }
    // [§68] T4 phase-coherent synthesis
    if (a == "--t4-phase-coherent") { t4_phase_coherent = true; continue; }
    if (a.rfind("--t4-alpha=", 0) == 0) { t4_alpha = std::stod(a.substr(11)); t4_phase_coherent = true; continue; }
    // [N1] Post-synthesis peak-level normalization.
    if (a == "--normalize-peak") {
      if (i + 1 >= argc) { std::cerr << "missing value for --normalize-peak\n"; usage(); return 1; }
      normalize_peak_dbfs    = std::stod(argv[++i]);
      normalize_peak_enabled = true;
      continue;
    }
    if (a.rfind("--normalize-peak=", 0) == 0) {
      normalize_peak_dbfs    = std::stod(a.substr(std::string("--normalize-peak=").size()));
      normalize_peak_enabled = true;
      continue;
    }
    if (a.rfind("--", 0) == 0) {
      std::cerr << "unknown option: " << a << "\n";
      usage();
      return 1;
    }
    pos_args.push_back(a);
  }

  // [V3 E3] --dbspec-refine-amp requires --dbspec-refine; auto-promote if used alone.
  // Documented: --dbspec-refine-amp without --dbspec-refine activates both.
  if (dbspec_refine_amp && !dbspec_refine) {
    dbspec_refine = true;
    std::cerr << "[V3 E3] info: --dbspec-refine-amp implies --dbspec-refine; auto-promoted.\n";
  }
  // [V3 E6] --grs-groups range check.
  if (grs_groups < 0) {
    std::cerr << "--grs-groups must be >= 1 (or omit to disable)\n";
    usage();
    return 1;
  }

  if (pos_args.size() > 10) {
    std::cerr << "too many positional arguments\n";
    usage();
    return 1;
  }

  const std::int64_t Nfft      = (pos_args.size() >= 1) ? std::stoll(pos_args[0]) : 1024;
  const std::int64_t hop       = (pos_args.size() >= 2) ? std::stoll(pos_args[1]) : 256;
  std::int64_t maxFrames       = (pos_args.size() >= 3) ? std::stoll(pos_args[2]) : 50;
  const double thresh_db       = (pos_args.size() >= 4) ? std::stod(pos_args[3])  : -300.0;

  if (frames_override > 0) {
    maxFrames = frames_override;
  }

  if ((synth_enabled || flt_enabled) && synth_out_wav.empty()) {
    std::cerr << "--synth/--flt requires --out-wav <path>\n";
    usage();
    return 1;
  }

  if (synth_skip < 0) {
    std::cerr << "--synth-skip must be >= 0\n";
    usage();
    return 1;
  }
  if (flt_enabled && flt_ir_path.empty()) {
    std::cerr << "--flt requires --flt-ir <fir.txt>\n";
    usage();
    return 1;
  }
  if (ndec <= 0) {
    std::cerr << "--ndec must be >= 1\n";
    usage();
    return 1;
  }
  // --maxoscs validation: reject zero/negative and pathologically small values.
  if (cli_maxoscs == 0 || (cli_maxoscs > 0 && cli_maxoscs < 10)) {
    std::cerr << "--maxoscs=" << cli_maxoscs
              << " invalid: must be >= 10 (or omit for default 300)\n";
    usage();
    return 1;
  }
  const std::int64_t maxoscs_val = (cli_maxoscs > 0) ? cli_maxoscs : 300;

  std::vector<double> x;
  int fs = 0;
  if (!parshl::ReadMonoAsDouble(inPath, x, fs)) return 2;

  // R1 OutPartials writer — open PREFIX.amp.wav / PREFIX.frq.wav if --out-partials was given.
  // Frame rate = round(Fs / hop).  Channel count = maxoscs_val (oscillator bank size).
  // Opened here so that any file-creation errors are reported before the analysis loop starts.
  parshl::PartialsWriter partials_writer;
  if (!out_partials_prefix.empty()) {
    const int frame_rate_hz = static_cast<int>(std::round(static_cast<double>(fs) / static_cast<double>(hop)));
    if (!partials_writer.open(out_partials_prefix, static_cast<int>(maxoscs_val), frame_rate_hz)) {
      std::cerr << "[OutPartials] aborting due to file-open error\n";
      return 5;
    }
  }

  const double binWidthHz = (Nfft > 0) ? (static_cast<double>(fs) / static_cast<double>(Nfft)) : 0.0;

  const double MinSepHz     = (pos_args.size() >= 5) ? std::stod(pos_args[4]) : binWidthHz;
  const double Hyst         = (pos_args.size() >= 6) ? std::stod(pos_args[5]) : 0.0;
  const std::int64_t MinWid = (pos_args.size() >= 7) ? std::stoll(pos_args[6]) : 1;
  const bool tracePartials  = (pos_args.size() >= 8) ? (std::stoll(pos_args[7]) != 0) : false;

  // Backward-compatible parsing:
  // - If only one extra arg after tracePartials: it's Nx
  // - If two extra args: [WinType] [Nx]
  int wintype = parshl::Hamming;
  std::int64_t Nx = Nfft;

  if (pos_args.size() == 9) {
    Nx = std::stoll(pos_args[8]);
  } else if (pos_args.size() == 10) {
    wintype = (int)std::stoll(pos_args[8]);
    Nx = std::stoll(pos_args[9]);
  }

  std::cout << "Read " << x.size() << " samples, fs=" << fs << "\n";
  std::cout << "Nfft=" << Nfft << " Nx=" << Nx << " hop=" << hop << " frames=" << maxFrames << "\n";
  std::cout << std::fixed << std::setprecision(6)
            << "binWidthHz=" << binWidthHz
            << " MinSepHz=" << MinSepHz
            << " Hyst=" << Hyst
            << " MinWid=" << MinWid
            << " thresh_db=" << thresh_db
            << " tracePartials=" << (tracePartials ? 1 : 0)
            << " DFmax1=" << (cli_dfmax1 >= 0.0 ? cli_dfmax1 : 20.0)
            << " DFmax2=" << (cli_dfmax2 >= 0.0 ? cli_dfmax2 : 200.0)
            << "\n";

  // STFT runner (now supports Nx != Nfft, with zero-padding)
  parshl::StftFrameRunner stft((int)Nfft, (int)hop, wintype, (int)Nx);

  // [V3 EXPERIMENT — Origin Type C]
  // Experiment:   SigScl analytic (E1)
  // Flag:         --sigscl-analytic
  // Origin:       Parshl-source.txt ~L944: "IS SIGSCL CORRECT? ? (Synth scaling)"
  // Source B:     DFT real cosine: |X[k]| = (A/2)·Σw  →  SigScl = 2 / sum(w[n])
  // Summary:      Replace SAIL heuristic 4/Nx with correct analytic formula.
  //               For Hann window: 4/Nfft = 2/Σw exactly (SAIL accidentally correct).
  //               For other windows (Hamming, Dolph-Chebyshev, …) factor 2 still holds.
  //               Only affects synthesis amplitude; tracker unchanged.
  // Contract:     Experimental only.  Disabled by default.
  //               V2 baseline bit-identical when flag is off.
  const double sail_sigscl = 4.0 / static_cast<double>(stft.Nx); // SAIL: SigScl <- 4/Nx
  double sum_window = 0.0;
  for (double w : stft.win) sum_window += w;
  const double analytic_sigscl = (sum_window > 0.0) ? 2.0 / sum_window : sail_sigscl;
  const double SigScl = sigscl_analytic ? analytic_sigscl : sail_sigscl;
  // end [V3 EXPERIMENT E1 --sigscl-analytic]
  parshl::FftwC2R ifft((int)Nfft);
  parshl::FftwR2C flt_fft((int)Nfft);
  const int nfft_i = static_cast<int>(Nfft);

  const std::int64_t Nspec = (Nfft/2 + 1);   // SAIL L307: Nspec <- (Nfft DIV 2)+1

  std::vector<double> flt_Hr;
  std::vector<double> flt_Hi;
  int flt_Nh = 0;
  if (flt_enabled) {
    std::vector<double> flt_ic;
    if (!load_fir_coeffs(flt_ir_path, flt_ic)) {
      std::cerr << "failed loading FIR from --flt-ir: " << flt_ir_path << "\n";
      return 1;
    }
    flt_Nh = static_cast<int>(flt_ic.size()); // SAIL: Nh <- Ni
    for (int i = 0; i < nfft_i; ++i) flt_fft.in[static_cast<std::size_t>(i)] = 0.0;
    for (int i = 0; i < flt_Nh && i < nfft_i; ++i) {
      flt_fft.in[static_cast<std::size_t>(i)] = flt_ic[static_cast<std::size_t>(i)];
    }
    flt_fft.exec();
    flt_Hr.resize(static_cast<std::size_t>(Nspec));
    flt_Hi.resize(static_cast<std::size_t>(Nspec));
    for (std::int64_t i = 0; i < Nspec; ++i) {
      flt_Hr[static_cast<std::size_t>(i)] = flt_fft.out[static_cast<std::size_t>(i)][0];
      flt_Hi[static_cast<std::size_t>(i)] = flt_fft.out[static_cast<std::size_t>(i)][1];
    }
  }

  // Tracker params (only for observing UpdateMap)
  parshl::TrackerParams tp;
  // Bank size: controlled by --maxoscs (default 300). MaxLins matches MaxOscs
  // (both govern the same oscillator pool; keep them in sync).
  tp.MaxOscs = maxoscs_val;
  tp.MaxLins = maxoscs_val;
  tp.Fs = static_cast<double>(fs);
  tp.Fc1 = 0.0;
  tp.Fc2 = fs * 0.5;

  // Placeholder values; in real PARSHL they come from α-D/α-U (DFmax1/DFmax2).
  // Driver defaults are kept intact when --dfmax1/--dfmax2 are not passed.
  // These flags only avoid hardcoding in the driver for tracker exploration;
  // they do not modify the UpdateMap/get_closest_frq/dfmax algorithm.
  tp.DFmax1 = (cli_dfmax1 >= 0.0) ? cli_dfmax1 : 20.0;
  tp.DFmax2 = (cli_dfmax2 >= 0.0) ? cli_dfmax2 : 200.0;
  tp.InstantRise = false;
  tp.UDtrace = false;
  tp.legacy_mode = legacy_mode;  // V2: recycling disabled when --legacy passed
  tp.dfmax_velocity_scale = legacy_mode ? 0.0 : cli_dfmax_vscale;  // V2 D7: adaptive DFmax (off in legacy mode)
  // D6-observe: disabled when --legacy; otherwise use CLI value (0 = off).
  tp.DAmax_observe = legacy_mode ? 0.0 : cli_damax_observe;
  // D6-gate: disabled when --legacy; otherwise use CLI value (0 = off).
  tp.DAmax_gate    = legacy_mode ? 0.0 : cli_damax_gate;
  // [V3 E2] Closest-free-oscillator dispatch: disabled when --legacy; otherwise use CLI value.
  tp.closest_osc_free = !legacy_mode && closest_osc_free;

  // D4: --reverse-analysis is a V2 feature and is suppressed in legacy mode.
  if (reverse_analysis && legacy_mode) {
    std::cerr << "[WARN] --reverse-analysis is ignored in --legacy mode (legacy preserves exact SAIL behaviour)\n";
    reverse_analysis = false;
  }

  parshl::ParshlTracker tracker(tp);

  // ---- Per-run statistics accumulators ----------------------------------------
  // Lifecycle (modes A/B: filled during main loop; mode C: filled from rev_frames).
  TrackAccum taccum{};
  // Birth-frame per oscillator slot for lifetime tracking (mode A/B main loop).
  std::vector<std::int64_t> osc_birth(static_cast<std::size_t>(maxoscs_val + 1), -1);
  // [§32B] DBspec refinement: accumulated |refined_frq - frq_before| per osc per frame.
  // Empty when --dbspec-refine is off (no allocation overhead).
  std::vector<double> dbspec_deltas;
  // [V3 E3] DBspec amplitude refinement: accumulated |amp_refined - amp_before| per osc per frame.
  // Empty when --dbspec-refine-amp is off.
  std::vector<double> dbspec_amp_deltas;
  // dfmax_exceeded from tracker_rev (mode C only; captured inside the D4 block).
  std::int64_t rev_dfmax_exceeded = 0;

  // ---- D4: Reverse-analysis pass (Serra & Smith 1990 §5 / SAIL L1389-1421) ----
  // Strategy from the paper:
  //   "The input sound is reversed in time to allow postponement of the attack
  //    analysis until the end. Working backwards through the sound, we dispatch
  //    a new oscillator on each new sinusoidal line which appears."
  // Implementation:
  //   1. Build a reversed copy of the input audio.
  //   2. Run a complete analysis pass (STFT → FindPartials → UpdateMap) on the
  //      ENTIRE reversed audio (sliding-window ring buffer of maxFrames entries).
  //      Walking from the END of the original (reversed pos 0) to the START,
  //      the ring buffer retains only the LAST maxFrames reversed frames, which
  //      correspond to the first maxFrames positions in FORWARD time (= attack).
  //   3. Reverse rev_frames[] so frame indices correspond to forward time.
  //   4. The main synthesis loop below uses those pre-computed parameters
  //      instead of a live tracker, ensuring attack transients receive
  //      oscillator parameters derived from their stable sustain (which the
  //      tracker "saw first" when walking the reversed audio).
  //
  // NOTE: earlier (buggy) version stopped the reversed pass after maxFrames
  // iterations, covering only the last maxFrames*hop samples of the original.
  // For signals longer than maxFrames*hop this missed the attack entirely.
  // The ring-buffer fix ensures the full audio is analysed regardless of length.
  //
  // When reverse_analysis == false, rev_frames is empty and the live tracker
  // path is used exactly as before — no change to existing behaviour.

  struct RevFrame {
    std::int64_t Noscs = 0;
    std::vector<double> OscAmp;   // 1-based, size MaxOscs+1
    std::vector<double> OscFrq;
    std::vector<double> PrvOscAmp;
    std::vector<double> PrvOscFrq;
    std::vector<std::int64_t> LinOfOsc;
    // Per-partial data (for phase init in the forward pass)
    std::int64_t Nlins = 0;
    std::vector<double> LinAmp_f;   // 1-based
    std::vector<double> LinFrq_f;
    std::vector<double> LinPhase_f; // spectral phase extracted during reverse pass
  };

  std::vector<RevFrame> rev_frames;

  if (reverse_analysis) {
    std::cout << "[D4] Running reverse-analysis pass on " << maxFrames << " frames...\n";

    // Build reversed audio (sample-level reversal).
    std::vector<double> x_rev(x.rbegin(), x.rend());

    parshl::StftFrameRunner stft_rev((int)Nfft, (int)hop, wintype, (int)Nx);
    parshl::ParshlTracker tracker_rev(tp);

    sail::Array1D<double> Xdb_rev(1, Nspec);
    sail::Array1D<double> LinAmpDB_rev(1, tp.MaxLins);
    sail::Array1D<double> LinFrq_rev_arr(1, tp.MaxLins);

    // Ring buffer of size maxFrames: keeps only the LAST maxFrames reversed frames,
    // which (after flipping) correspond to the first maxFrames positions in forward
    // time (= attack + early sustain). Older frames (= tail of original) are discarded
    // as they slide out the front of the deque.
    const std::int64_t total_rev_est =
      (static_cast<std::int64_t>(x_rev.size()) + hop - 1) / hop;
    std::cout << "[D4] Running reverse-analysis pass: est " << total_rev_est
              << " frames total, keeping last " << maxFrames << " for synthesis...\n";

    std::deque<RevFrame> rev_ring;
    std::int64_t pos_rev = 0;
    std::int64_t f_rev = 0;

    for (;;) {
      if (pos_rev >= (std::int64_t)x_rev.size()) break;

      stft_rev.compute_db(x_rev, (long)pos_rev, Xdb_rev);

      const std::int64_t Nlins_rev = parshl::FindPartials(
        Xdb_rev, LinAmpDB_rev, LinFrq_rev_arr,
        (double)fs, (std::int64_t)Nfft,
        tp.Fc1, tp.Fc2, MinSepHz, thresh_db, Hyst,
        tp.MaxLins, MinWid, /*trace=*/false
      );

      // dB → linear amplitude (SAIL L1728)
      const std::size_t linsz = static_cast<std::size_t>(tp.MaxLins + 1);
      std::vector<double> LA_rev(linsz, 0.0);
      std::vector<double> LF_rev(linsz, 0.0);
      std::vector<double> LP_rev(linsz, 0.0); // spectral phases
      for (std::int64_t i = 1; i <= Nlins_rev; ++i) {
        LF_rev[static_cast<std::size_t>(i)] = LinFrq_rev_arr[i];
        LA_rev[static_cast<std::size_t>(i)] = std::pow(10.0, LinAmpDB_rev[i] / 20.0) * SigScl;
      }

      // Complex peak interpolation (D1/D2/D3) — same V2 path as forward pass
      for (std::int64_t i = 1; i <= Nlins_rev; ++i) {
        const double bin_frac = LF_rev[static_cast<std::size_t>(i)]
                                * static_cast<double>(Nfft) / static_cast<double>(fs) + 1.0;
        const std::int64_t k = static_cast<std::int64_t>(std::llround(bin_frac));
        const double p = bin_frac - static_cast<double>(k);
        double re = 0.0, im = 0.0;
        parshl::InterpComplexAtPeak(stft_rev.lastS, Nspec, k, p, re, im);
        const double amp_interp = std::sqrt(re * re + im * im);
        if (amp_interp > 0.0) {
          LA_rev[static_cast<std::size_t>(i)] = amp_interp * SigScl;
          LP_rev[static_cast<std::size_t>(i)] = std::atan2(im, re);
        }
      }

      tracker_rev.update_map((long)(f_rev + 1), (long)Nlins_rev, LA_rev, LF_rev);

      // Snapshot this frame's synthesiser parameters.
      const auto& st = tracker_rev.state();
      RevFrame rf;
      rf.Noscs = st.Noscs;
      rf.OscAmp    = st.OscAmp;
      rf.OscFrq    = st.OscFrq;
      rf.PrvOscAmp = st.PrvOscAmp;
      rf.PrvOscFrq = st.PrvOscFrq;
      rf.LinOfOsc  = std::vector<std::int64_t>(st.LinOfOsc.begin(), st.LinOfOsc.end());
      rf.Nlins     = Nlins_rev;
      rf.LinAmp_f  = LA_rev;
      rf.LinFrq_f  = LF_rev;
      rf.LinPhase_f = LP_rev;

      // Ring buffer: evict oldest if full
      if (static_cast<std::int64_t>(rev_ring.size()) >= maxFrames) {
        rev_ring.pop_front();
      }
      rev_ring.push_back(std::move(rf));

      pos_rev += hop;
      ++f_rev;
    }

    // Transfer ring buffer → rev_frames (ring holds last maxFrames reversed frames
    // = frames nearest the START of the original audio = attack/early-sustain region)
    rev_frames.assign(rev_ring.begin(), rev_ring.end());

    // Flip so rev_frames[0] = last reversed frame = first forward frame.
    std::reverse(rev_frames.begin(), rev_frames.end());

    // Fix up Prv* fields: after reversing, frame[f].Prv* should equal frame[f-1].Cur*.
    // The OscAmp/OscFrq at frame[f] in forward time = OscAmp/OscFrq at rev_frames[f]
    // (already the current-frame values from the reversed tracker).
    // PrvOscAmp/PrvOscFrq for frame[f] = OscAmp/OscFrq from frame[f-1].
    // Re-derive them from the already-flipped sequence.
    for (std::size_t f = rev_frames.size(); f-- > 0; ) {
      if (f == 0) {
        // First forward frame: Prv = 0 (oscillators just born)
        std::fill(rev_frames[0].PrvOscAmp.begin(), rev_frames[0].PrvOscAmp.end(), 0.0);
        // Keep PrvOscFrq = OscFrq (OutAF rule: PrvOscAmp==0 → PrvOscFrq=OscFrq)
        rev_frames[0].PrvOscFrq = rev_frames[0].OscFrq;
      } else {
        // frame[f].Prv* = frame[f-1].Cur*
        rev_frames[f].PrvOscAmp = rev_frames[f - 1].OscAmp;
        rev_frames[f].PrvOscFrq = rev_frames[f - 1].OscFrq;
      }
    }

    std::cout << "[D4] Reverse pass complete: " << f_rev << " frames analysed, "
              << rev_frames.size() << " retained for synthesis.\n";

    // ---- Compute lifecycle stats from rev_frames (mode C synthesis source) -------
    // Build the parallel vectors needed by compute_revframes_lifecycle().
    {
      std::vector<std::pair<std::int64_t, std::vector<std::int64_t>>> frames_lof;
      std::vector<std::int64_t> noscs_v, nlins_v;
      frames_lof.reserve(rev_frames.size());
      noscs_v.reserve(rev_frames.size());
      nlins_v.reserve(rev_frames.size());
      for (const auto& rf : rev_frames) {
        frames_lof.emplace_back(rf.Noscs, rf.LinOfOsc);
        noscs_v.push_back(rf.Noscs);
        nlins_v.push_back(rf.Nlins);
      }
      taccum = compute_revframes_lifecycle(frames_lof, maxoscs_val, noscs_v, nlins_v);
    }
    // dfmax_exceeded from the reversed tracker (reflects frequency-jump stops
    // across the FULL reversed pass — covers attack-to-sustain transition).
    // lifecycle add/stop/recycle come from compute_revframes_lifecycle() which
    // operates on the RETAINED frames only (consistent with n_frames, Noscs_avg, etc.).
    rev_dfmax_exceeded   = tracker_rev.stats().dfmax_exceeded;
    // Note: taccum.add_count/stop_count/recycle_count are already set by
    // compute_revframes_lifecycle() above; do NOT overwrite with the full-pass counts.
  }
  // ---- end D4 reverse pass ----
  sail::Array1D<double> Xdb(1, Nspec);

  sail::Array1D<double> Peaks(1, 600);
  sail::Array1D<double> Locs (1, 600);

  sail::Array1D<double> LinAmpDB(1, tp.MaxLins);
  sail::Array1D<double> LinFrq(1, tp.MaxLins);

  parshl::ParshlSynthState synth_state;
  std::vector<float> OutBuf;
  std::vector<float> out_buf;
  std::size_t out_written = 0;
  int Bp = 0;
  int decim_P = 1;
  bool decim_first_call = true;
  const bool synth_log = synth_log_enabled();
  const int synth_nfft = static_cast<int>(Nfft);
  const int Ndec = ndec;
  const bool out_enabled = (synth_enabled || flt_enabled);
  if (out_enabled) {
    OutBuf.assign(static_cast<std::size_t>(2 * synth_nfft), 0.0f);
    out_buf.clear();
    out_buf.reserve(static_cast<std::size_t>(std::max<std::int64_t>(0, maxFrames) * hop));
  }

  std::int64_t pos = 0;

  // SAIL main analysis loop "RI" (Parshl-source.txt L1690–L1780)
  // Per frame: STFT → FindPartials → LinAmp conversion → UpdateMap → Synthesize → OLA flush
  for (std::int64_t f = 0; f < maxFrames; ++f) {
    // SAIL L1695–L1717: ARRCLR + read + window + !FFA + XmagDB
    stft.compute_db(x, (long)pos, Xdb);

    // ---- FindPartials (SAIL L1726–L1727) ----
    const double Fc1 = tp.Fc1;
    const double Fc2 = tp.Fc2;

    const std::int64_t NpartialsReq = tp.MaxLins;

    const std::int64_t Nlins = parshl::FindPartials(
      Xdb,
      LinAmpDB,
      LinFrq,
      (double)fs,
      (std::int64_t)Nfft,
      Fc1,
      Fc2,
      MinSepHz,
      thresh_db,
      Hyst,
      NpartialsReq,
      MinWid,
      /*trace=*/tracePartials
    );

    const std::int64_t I1 = 1;
    const std::int64_t I2 = Nspec;

    std::cout << "FindPartials: nfound=" << Nlins
              << " I1=" << I1 << " I2=" << I2;

    if (tracePartials) {
      const std::int64_t BinInt = BinInt_from_MinSepHz(MinSepHz, (double)fs, Nfft);
      std::cout << " BinInt=" << BinInt;
    }

    std::cout << " (Fc1=" << Fc1
              << " Fc2=" << Fc2
              << " MinSep=" << MinSepHz
              << " Hyst=" << Hyst
              << " thresh_db=" << thresh_db
              << ")\n";

    if (tracePartials) {
      const std::int64_t nprint = std::min<std::int64_t>(Nlins, 10);
      for (std::int64_t i = 1; i <= nprint; ++i) {
        std::cout << "  lin " << i
                  << ": ampDB=" << LinAmpDB[i]
                  << " frqHz=" << LinFrq[i] << "\n";
      }
    }

    // SAIL L1728: LinAmp[i] <- 10^(LinAmpDB[i]/20) * SigScl   (dB → linear amplitude)
    // LinFrq_v: vector copy of LinFrq (sail::Array1D) needed by tracker.update_map (std::vector).
    std::vector<double> LinAmp(static_cast<std::size_t>(tp.MaxLins + 1), 0.0);
    std::vector<double> LinFrq_v(static_cast<std::size_t>(tp.MaxLins + 1), 0.0);
    // V2: spectral phase per partial (from complex interpolation). Zero in legacy mode.
    std::vector<double> LinPhase(static_cast<std::size_t>(tp.MaxLins + 1), 0.0);
    for (std::int64_t i = 1; i <= Nlins; ++i) {
      LinFrq_v[i] = LinFrq[i];
      LinAmp[i] = std::pow(10.0, LinAmpDB[i] / 20.0) * SigScl;
    }

    // --- V2: Step 4 — Complex Peak Interpolation (Serra & Smith 1990 §5) ----------
    // Replace integer-bin amplitude with parabolic interpolation of Re/Im at the
    // fractional peak location. Also extracts spectral phase for new oscillators.
    // --legacy flag: skip this block → exact V1 (SAIL) behaviour.
    if (!legacy_mode) {
      for (std::int64_t i = 1; i <= Nlins; ++i) {
        // Recover fractional bin from Hz (inverse of FindPartials bin→Hz, L470):
        //   LinFrq[i] = Fs*(bin_frac - 1)/Nfft  →  bin_frac = i*Fs/Nfft + 1
        const double bin_frac = LinFrq_v[i] * static_cast<double>(Nfft)
                                / static_cast<double>(fs) + 1.0;
        const std::int64_t k = static_cast<std::int64_t>(std::llround(bin_frac));
        const double p = bin_frac - static_cast<double>(k);
        double re = 0.0, im = 0.0;
        parshl::InterpComplexAtPeak(stft.lastS, Nspec, k, p, re, im);
        const double amp_interp = std::sqrt(re * re + im * im);
        if (amp_interp > 0.0) {
          // Replace dB-at-integer-bin amplitude with complex-interpolated magnitude.
          LinAmp[i]   = amp_interp * SigScl;
          // Save spectral phase (radians, [-π,π]) for oscillator init below.
          LinPhase[static_cast<std::size_t>(i)] = std::atan2(im, re);
        }
        // amp_interp == 0 is a degenerate edge case; keep dB-based value.
      }
    }
    // --------------------------------------------------------------------------

    // [§65] freq_error: |f_interpolated - f_integer_bin| per partial.
    // = |Qinterp_offset| * Fs/Nfft (Hz).  Measures the sub-bin correction applied
    // by FindPeaks parabolic interpolation; large values = peaks far from bin centres.
    for (std::int64_t i = 1; i <= Nlins; ++i) {
      const double bin_frac = LinFrq_v[static_cast<std::size_t>(i)]
                              * static_cast<double>(Nfft) / static_cast<double>(fs) + 1.0;
      const double fe = std::abs(bin_frac - std::round(bin_frac))
                        * static_cast<double>(fs) / static_cast<double>(Nfft);
      diag_freq_err_sum += fe;
      diag_freq_err_max  = std::max(diag_freq_err_max, fe);
      ++diag_freq_err_n;
    }

    if (trk_run_log_enabled()) {
      const auto& st_before = tracker.state();
      const std::int64_t want = std::min<std::int64_t>(Nlins, tp.MaxOscs);
      const std::int64_t Nadd = want - st_before.Noscs;
      std::cout << "[TRK_RUN] frame=" << f
                << " frame1=" << (f + 1)
                << " Nlins=" << Nlins
                << " Noscs_before=" << st_before.Noscs
                << " want=" << want
                << " Nadd=" << Nadd
                << " InstantRise=" << (tp.InstantRise ? 1 : 0)
                << " fs=" << fs
                << " Nx=" << Nx
                << " hop=" << hop
                << " Nfft=" << Nfft
                << " binWidthHz=" << binWidthHz
                << " MinSepHz=" << MinSepHz
                << " Hyst=" << Hyst
                << " thresh_db=" << thresh_db
                << " DFmax1=" << tp.DFmax1
                << " DFmax2=" << tp.DFmax2
                << " Fc1=" << tp.Fc1
                << " Fc2=" << tp.Fc2
                << "\n";
    }

    // SAIL UpdateMap call (tracker)
    // D4: when reverse_analysis is active, the live tracker is still run here
    // (so that console output / track stats are consistent), but synthesis will
    // use rev_frames[f] parameters instead of the live tracker state.
    // [V3 EXPERIMENT — Origin Type B]
    // Experiment:   GRS spectral group reduction (E6)
    // Flag:         --grs-groups=N
    // Origin B:     GRS comment (Parshl-source.txt L1432–L1436):
    //               "GRS suggested doing groups reduction on the mag² spectrum,
    //                letting each group behave as a 'virtual line'."
    // Algorithm:    Divide [Fc1,Fc2] into N equal-bandwidth groups.
    //               Per group: sum mag² of all STFT bins in the group;
    //               centroid_frq = sum(frq_b * mag2_b) / sum(mag2_b);
    //               virtual_amp   = sqrt(sum_mag2) * SigScl.
    //               Skips groups with sum_mag2 == 0. Virtual lines overwrite
    //               LinFrq_v[1..Nlins_virt] / LinAmp[1..Nlins_virt].
    // Metric:       grs_group_collision_count: cumulative (frame×group) pairs
    //               where ≥2 FindPartials peaks fall in the same group.
    // Note:         FindPartials still runs; Nlins used for taccum.total_nlins.
    //               LinPhase remains indexed 1..Nlins; virtual lines get phase=0
    //               (SAIL default for new oscillators). Synthesis unaffected when off.
    // Contract:     Experimental. Disabled by default. V2 baseline unchanged when off.
    std::int64_t Nlins_eff = Nlins;
    if (grs_groups > 0 && !reverse_analysis) {
      const double Fc2_grs = (tp.Fc2 > 0.0) ? tp.Fc2 : (static_cast<double>(fs) * 0.5);
      const double group_bw = (Fc2_grs - Fc1) / static_cast<double>(grs_groups);
      // ---- Collision counting: FindPartials peaks per group. -------------------
      std::vector<std::int64_t> group_peak_count(
          static_cast<std::size_t>(grs_groups), 0LL);
      if (group_bw > 0.0) {
        for (std::int64_t i = 1; i <= Nlins; ++i) {
          const double f_hz = LinFrq_v[static_cast<std::size_t>(i)];
          if (f_hz >= Fc1 && f_hz < Fc2_grs) {
            const std::int64_t gc = std::max<std::int64_t>(0,
                std::min(grs_groups - 1,
                    static_cast<std::int64_t>(std::floor((f_hz - Fc1) / group_bw))));
            ++group_peak_count[static_cast<std::size_t>(gc)];
          }
        }
      }
      for (std::int64_t g = 0; g < grs_groups; ++g) {
        if (group_peak_count[static_cast<std::size_t>(g)] >= 2)
          ++grs_group_collision_count;
      }
      // ---- Build virtual lines from Xdb mag² spectrum. -----------------------
      // Extend LinFrq_v / LinAmp if grs_groups > current allocation.
      const std::size_t needed = static_cast<std::size_t>(grs_groups + 1);
      if (LinFrq_v.size() < needed) LinFrq_v.assign(needed, 0.0);
      if (LinAmp.size()  < needed) LinAmp.assign(needed, 0.0);
      std::int64_t Nlins_virt = 0;
      for (std::int64_t g = 0; g < grs_groups; ++g) {
        const double f_lo = Fc1 + static_cast<double>(g)     * group_bw;
        const double f_hi = Fc1 + static_cast<double>(g + 1) * group_bw;
        // 1-based bin range: bin b → frequency (b−1)·fs/Nfft (Parshl-source.txt L470).
        const long b_lo = std::max(1L,
            static_cast<long>(std::floor(f_lo * static_cast<double>(Nfft)
                                         / static_cast<double>(fs))) + 1);
        const long b_hi = std::min((long)Nspec,
            static_cast<long>(std::floor(f_hi * static_cast<double>(Nfft)
                                         / static_cast<double>(fs))));
        if (b_lo > b_hi) continue;  // empty bin range at spectral boundaries
        double sum_mag2 = 0.0;
        double sum_wfrq = 0.0;
        for (long b = b_lo; b <= b_hi; ++b) {
          // Xdb[b] is in dB (DBscl space): mag² = 10^(db/10).
          const double mag2  = std::pow(10.0, Xdb[b] / 10.0);
          const double frq_b = static_cast<double>(fs)
              * (static_cast<double>(b) - 1.0) / static_cast<double>(Nfft);
          sum_mag2 += mag2;
          sum_wfrq += mag2 * frq_b;
        }
        if (sum_mag2 <= 0.0) continue;  // silent group — skip
        ++Nlins_virt;
        LinFrq_v[static_cast<std::size_t>(Nlins_virt)] = sum_wfrq / sum_mag2;
        LinAmp  [static_cast<std::size_t>(Nlins_virt)] = std::sqrt(sum_mag2) * SigScl;
      }
      Nlins_eff = Nlins_virt;
    }
    // end [V3 EXPERIMENT E6 --grs-groups]

    // [V4 T1 / V3 V2 dispatch]
    // When --t1-peak-claim is active (and not reverse_analysis), run the
    // peak-claiming tracker and inject its result back into tracker.state()
    // so that all downstream synthesis and lifecycle-counting code is
    // unchanged.  Otherwise fall back to the PARSHL V2 greedy tracker.
    if (t1_peak_claim && !reverse_analysis) {
      parshl::TrackerState t1_new_state;
      const auto t1s = parshl::t1_update(
          tracker.state(), tp, (long)(f + 1),
          (long)Nlins_eff, LinAmp, LinFrq_v,
          t1_squelch_hold, t1_hold_counts,
          t1_amp_hysteresis,  // [§64]
          t1_new_state);
      t1_stats_accum.t1_created             += t1s.t1_created;
      t1_stats_accum.t1_assigned            += t1s.t1_assigned;
      t1_stats_accum.t1_terminated          += t1s.t1_terminated;
      t1_stats_accum.t1_hold_survivals      += t1s.t1_hold_survivals;
      t1_stats_accum.t1_hold_terminations   += t1s.t1_hold_terminations;
      tracker.overwrite_state(std::move(t1_new_state));
    } else if (t2_global_assign && !reverse_analysis) {
      // [V4 T2] Global assignment tracker: min-cost bipartite matching per frame.
      parshl::TrackerState t2_new_state;
      const auto t2s = parshl::t2_update(
          tracker.state(), tp, (long)(f + 1),
          (long)Nlins_eff, LinAmp, LinFrq_v,
          t2_trajectory_cost ? t2_wid : 0.0,   // [§60] pass wid (0 = disabled)
          t2_rank_cost ? t2_wrank : 0.0,        // [§61] pass wrank (0 = disabled)
          t2_last_peak_indices,                 // [§60] per-slot trajectory state
          t2_phase_cost ? t2_wphi : 0.0,        // [§66B] phase weight (0 = disabled)
          static_cast<double>(hop) / static_cast<double>(fs), // [§66B] dt = hop/fs
          LinPhase,                             // [§66B] spectral phases, 1-based
          t2_prv_phase,                         // [§66B] per-slot previous phase
          t2_new_state);
      t2_stats_accum.assignments       += t2s.assignments;
      t2_stats_accum.creations         += t2s.creations;
      t2_stats_accum.terminations      += t2s.terminations;
      t2_stats_accum.identity_switches += t2s.identity_switches;  // [§60]
      t2_stats_accum.phase_eligible    += t2s.phase_eligible;     // [§66B]
      t2_stats_accum.phase_cost_sum    += t2s.phase_cost_sum;     // [§66B]
      tracker.overwrite_state(std::move(t2_new_state));
    } else {
      tracker.update_map((long)(f + 1), (long)Nlins_eff, LinAmp, LinFrq_v);
    }

    // [V3 EXPERIMENT — Origin Type B+C]
    // Experiment:   DBspec matching (E5)
    // Flag:         --dbspec-match
    // Origin B:     JOS 8-JUL-85 (Parshl-source.txt L1425–L1432):
    //               "Consider allowing Oscs to find their peaks in DBspec instead of LinFrqs."
    // Origin C:     Conservative Variant A: maintain FindPartials + LinFrq + update_map.
    //               After update_map() assigns LinFrq candidates, search Xdb within
    //               ±DFmax_bins around the PREDICTED frequency (PrvOscFrq = last-frame OscFrq).
    //               Replace the LinFrq match with the DBspec peak only if it is strictly
    //               closer to the predicted frequency (|dbspec_frq - pred| < |linfrq - pred|).
    // Interaction:  Runs BEFORE --dbspec-refine (coarse correction → fine refinement).
    //               New/recycled oscillators are automatically excluded: OutAF rule sets
    //               PrvOscFrq = OscFrq when PrvOscAmp == 0, so dist_linfrq == 0 → no override.
    // Contract:     Experimental. Disabled by default. V2 baseline unchanged when off.
    if (dbspec_match && !reverse_analysis) {
      const double Fc2_hz = (tp.Fc2 > 0.0) ? tp.Fc2 : (static_cast<double>(fs) * 0.5);
      const std::int64_t lim = std::min(maxoscs_val,
          static_cast<std::int64_t>(tracker.state().LinOfOsc.size()) - 1);
      const auto& dm_st = tracker.state();
      for (std::int64_t osc = 1; osc <= lim; ++osc) {
        if (!tracker.osc_is_on(osc)) continue;
        const std::size_t oi = static_cast<std::size_t>(osc);
        // predicted_frq: frequency BEFORE this frame's match (previous frame's OscFrq).
        // For new/recycled oscillators OutAF sets PrvOscFrq = OscFrq → dist_linfrq == 0
        // → criterion never satisfied → they are silently excluded.
        const double predicted_frq = dm_st.PrvOscFrq[oi];
        const double matched_frq   = dm_st.OscFrq[oi];  // LinFrq candidate from update_map
        // Per-oscillator DFmax window: linear model DFmax(f) = DFmax1 + (DFmax2-DFmax1)*f/Fc2
        const double dfy = tp.DFmax1
            + (tp.DFmax2 - tp.DFmax1) * predicted_frq / Fc2_hz;
        const long dfmax_bins = std::max(1L, static_cast<long>(std::round(
            dfy * static_cast<double>(Nfft) / static_cast<double>(fs))));
        // bin_center: 1-based bin nearest predicted_frq (Parshl-source.txt L470 convention)
        const long bin_c = static_cast<long>(
            predicted_frq * static_cast<double>(Nfft) / static_cast<double>(fs) + 1.5);
        const long lo = std::max(1L,        bin_c - dfmax_bins);
        const long hi = std::min((long)Nspec, bin_c + dfmax_bins);
        // Local max search in Xdb (read-only — no clobber)
        long best_bin = (bin_c >= 1 && bin_c <= (long)Nspec) ? bin_c : lo;
        double best_val = Xdb[best_bin];
        for (long b = lo; b <= hi; ++b) {
          if (Xdb[b] > best_val) { best_val = Xdb[b]; best_bin = b; }
        }
        // Sub-bin interpolation via Qinterp on dB (SAIL L396 / parshl_findpeaks.cpp)
        const long Lb = std::max(1L,        best_bin - 1);
        const long Rb = std::min((long)Nspec, best_bin + 1);
        const double offset = parshl::Qinterp(Xdb[Lb], Xdb[best_bin], Xdb[Rb]);
        double dbspec_frq = static_cast<double>(fs)
            * (static_cast<double>(best_bin) + offset - 1.0)
            / static_cast<double>(Nfft);
        // Defensive clip to [Fc1, Fc2]
        dbspec_frq = std::max(Fc1, std::min(Fc2_hz, dbspec_frq));
        // Criterion: replace only if DBspec peak is STRICTLY closer to predicted_frq
        const double dist_dbspec = std::fabs(dbspec_frq  - predicted_frq);
        const double dist_linfrq = std::fabs(matched_frq - predicted_frq);
        if (dist_dbspec < dist_linfrq) {
          tracker.set_osc_frq(osc, dbspec_frq);
          ++dbspec_match_count;
        }
      }
    }
    // end [V3 EXPERIMENT E5 --dbspec-match]

    // [§32B] --dbspec-refine: post-match local frequency refinement in DBspec.
    // JOS 8-JUL-85 (Parshl-source.txt L1425–L1432):
    //   "Consider allowing Oscs to find their peaks in DBspec instead of LinFrqs.
    //    Operation would be exactly as now, except the raw DB mag spectrum is searched
    //    for a local max nearest the current osc frq."
    // V2 conservative approximation: matching is unchanged; only OscFrq is refined.
    // Xdb is read-only (no clobber, no peak removal).
    // Not applied to D4/reverse-analysis in this session.
    if (dbspec_refine && !reverse_analysis) {
      // half_win_bins = max(1, round(MinSepHz/2 * Nfft/Fs))
      // Keeps the search window strictly local: at most half the partial separation.
      const long half_win = std::max(1L,
          static_cast<long>(std::round((MinSepHz / 2.0)
                                       * static_cast<double>(Nfft)
                                       / static_cast<double>(fs))));
      const std::int64_t lim = std::min(maxoscs_val,
          static_cast<std::int64_t>(tracker.state().LinOfOsc.size()) - 1);
      for (std::int64_t osc = 1; osc <= lim; ++osc) {
        if (!tracker.osc_is_on(osc)) continue;
        const double frq_before = tracker.osc_frq_hz(osc);
        // bin_c: 1-based bin index closest to frq_before (Parshl-source.txt L470 convention)
        const long bin_c = static_cast<long>(
            frq_before * static_cast<double>(Nfft) / static_cast<double>(fs) + 1.5);
        const long lo = std::max(1L,       bin_c - half_win);
        const long hi = std::min((long)Nspec, bin_c + half_win);
        // Local max search in Xdb (read-only, same DBscl space as FindPartials)
        long   best_bin = (bin_c >= 1 && bin_c <= (long)Nspec) ? bin_c : lo;
        double best_val = Xdb[best_bin];
        for (long b = lo; b <= hi; ++b) {
          if (Xdb[b] > best_val) { best_val = Xdb[b]; best_bin = b; }
        }
        // Sub-bin interpolation via Qinterp on dB (SAIL L396 / parshl_findpeaks.cpp)
        const long Lb = std::max(1L,       best_bin - 1);
        const long Rb = std::min((long)Nspec, best_bin + 1);
        const double offset = parshl::Qinterp(Xdb[Lb], Xdb[best_bin], Xdb[Rb]);
        double refined_frq = static_cast<double>(fs)
            * (static_cast<double>(best_bin) + offset - 1.0)
            / static_cast<double>(Nfft);
        // Defensive clip to [Fc1, Fc2] (same bounds as FindPartials)
        refined_frq = std::max(Fc1, std::min(Fc2, refined_frq));
        tracker.set_osc_frq(osc, refined_frq);
        dbspec_deltas.push_back(std::abs(refined_frq - frq_before));
        // [V3 EXPERIMENT — Origin Type B+C]
        // Experiment:   DBspec amplitude refinement (E3)
        // Flag:         --dbspec-refine-amp
        // Origin B:     JOS 8-JUL-85 (Parshl-source.txt L1425): "Consider allowing Oscs
        //               to find their peaks in DBspec instead of LinFrqs."
        // Origin C:     Analytical: Qinterp peak value → amp via 10^(db/20)*SigScl.
        // Summary:      NOT full DBspec matching.  Conservative post-match amplitude update
        //               using the same Xdb local peak already identified for freq refinement.
        //               db_peak = Y0 + offset*(Yp1-Ym1)/2  (parabolic value at interpolated bin)
        //               amp_refined = 10^(db_peak/20) * SigScl  (identical to LinAmp pipeline)
        // Contract:     Experimental. Disabled by default. V2 baseline unchanged when off.
        if (dbspec_refine_amp) {
          const double db_peak = Xdb[best_bin]
              + offset * (Xdb[Rb] - Xdb[Lb]) / 2.0;
          const double amp_refined = std::pow(10.0, db_peak / 20.0) * SigScl;
          const double amp_before  = tracker.osc_amp_lin(osc);
          tracker.set_osc_amp(osc, amp_refined);
          dbspec_amp_deltas.push_back(std::abs(amp_refined - amp_before));
        }
        // end [V3 EXPERIMENT E3 --dbspec-refine-amp]
      }
    }

    // ---- Per-frame lifecycle / stats accumulation --------------------------------
    // For mode C (reverse_analysis), taccum was pre-filled from rev_frames;
    // the live forward tracker is not used for synthesis so we skip lifecycle here.
    // For modes A/B, accumulate from the live tracker state after update_map().
    if (!reverse_analysis) {
      const auto& st_acc = tracker.state();
      ++taccum.n_frames;
      taccum.total_nlins += Nlins;
      taccum.max_nlins = std::max(taccum.max_nlins, Nlins);
      taccum.total_noscs += st_acc.Noscs;
      taccum.max_noscs = std::max(taccum.max_noscs, st_acc.Noscs);

      // Detect lifecycle transitions: PrvLinOfOsc = frame f-1 state,
      //                               LinOfOsc    = frame f state.
      const std::int64_t lim_acc = std::min(maxoscs_val,
                                    static_cast<std::int64_t>(st_acc.LinOfOsc.size()) - 1);
      for (std::int64_t osc = 1; osc <= lim_acc; ++osc) {
        const std::size_t oi = static_cast<std::size_t>(osc);
        if (oi >= st_acc.PrvLinOfOsc.size() || oi >= st_acc.LinOfOsc.size()) continue;
        const std::int64_t prv = st_acc.PrvLinOfOsc[oi];
        const std::int64_t cur = st_acc.LinOfOsc[oi];
        if (prv <= 0 && cur > 0) {     // trajectory start
          if (prv < 0) ++taccum.recycle_count;
          else         ++taccum.add_count;
          osc_birth[oi] = f;
        }
        if (prv > 0 && cur <= 0) {     // trajectory end
          ++taccum.stop_count;
          if (osc_birth[oi] >= 0) {
            taccum.lifetimes.push_back(f - osc_birth[oi]);
            osc_birth[oi] = -1;
          }
        }
      }
    }

    // [§65] amp_error: |LinAmp[lin] - OscAmp[osc]| for each matched oscillator.
    // Quantifies analysis-to-synthesis amplitude deviation (T1b hysteresis,
    // dbspec-refine-amp, etc.).  Zero for C1/C2 when no post-processing modifies OscAmp.
    if (!reverse_analysis) {
      const auto& st_ae = tracker.state();
      const std::int64_t lim_ae = std::min(maxoscs_val,
          static_cast<std::int64_t>(st_ae.LinOfOsc.size()) - 1);
      for (std::int64_t osc = 1; osc <= lim_ae; ++osc) {
        const std::size_t oi = static_cast<std::size_t>(osc);
        if (oi >= st_ae.LinOfOsc.size()) continue;
        const std::int64_t lin = st_ae.LinOfOsc[oi];
        if (lin < 1 || lin > Nlins_eff
            || static_cast<std::size_t>(lin) >= LinAmp.size()
            || oi >= st_ae.OscAmp.size()) continue;
        diag_amp_err_sum += std::abs(st_ae.OscAmp[oi]
                                     - LinAmp[static_cast<std::size_t>(lin)]);
        ++diag_amp_err_n;
      }
    }

    // --- V2: initialize OscPhs for newly-started oscillators from spectral phase ----
    // SAIL: OscPhs[osc] <- 0  for all new oscillators (Frame 1 init + Nadd adds).
    // V2:   OscPhs[osc] = theta / (2π) * SinSiz  where theta = atan2(Im, Re) at peak.
    // --legacy flag: skip → OscPhs stays 0 (SAIL/V1 behaviour).
    // D4: when reverse_analysis is active, use the LinPhase from the reversed frame
    //     (rev_frames[f]) for newly-started oscillators, because those are the
    //     spectral phases observed at the correct forward-time position.
    if (!legacy_mode && synth_enabled) {
      // Choose source of phase data: D4 reversed frame or current forward frame.
      const std::vector<double>* phase_src = &LinPhase;
      std::vector<double> rev_phase_local; // only used in D4 path
      const auto* src_st = &tracker.state(); // live tracker state (for LinOfOsc)

      if (reverse_analysis && static_cast<std::size_t>(f) < rev_frames.size()) {
        rev_phase_local = rev_frames[static_cast<std::size_t>(f)].LinPhase_f;
        phase_src = &rev_phase_local;
      }

      const auto& st = *src_st;
      // Ensure OscPhs is large enough before we write into it.
      const std::size_t phs_need = static_cast<std::size_t>(tp.MaxOscs + 1);
      if (synth_state.OscPhs.size() < phs_need) {
        synth_state.OscPhs.resize(phs_need, 0.0);
      }
      constexpr double two_pi = 2.0 * std::numbers::pi;
      const double SinSiz_d = static_cast<double>(parshl::ParshlSynthState::SinSiz);
      for (std::int64_t osc = 1; osc <= st.Noscs; ++osc) {
        const std::size_t osc_sz = static_cast<std::size_t>(osc);
        if (osc_sz >= st.LinOfOsc.size()) continue;
        // Newly-started oscillator: was Free (PrvLinOfOsc==0) before, is On (>0) now.
        // Frame 1 (f==0): all allocated oscillators are new (SAIL sets PrvLinOfOsc=LinOfOsc).
        const bool was_free = (osc_sz < st.PrvLinOfOsc.size())
                              && (st.PrvLinOfOsc[osc_sz] == 0);
        // V2 recycling: also init phase for SQUELCHED→ON oscillators that were fully silent.
        const bool was_squelched_silent =
            (osc_sz < st.PrvLinOfOsc.size())
            && (st.PrvLinOfOsc[osc_sz] < 0)
            && (osc_sz < st.PrvOscAmp.size())
            && (st.PrvOscAmp[osc_sz] == 0.0);
        const bool is_on    = st.LinOfOsc[osc_sz] > 0;
        const bool is_new   = (f == 0) ? is_on : ((was_free || was_squelched_silent) && is_on);
        if (!is_new) continue;
        const std::int64_t lin = st.LinOfOsc[osc_sz];
        if (lin < 1 || lin >= static_cast<std::int64_t>(phase_src->size())) continue;
        const double theta = (*phase_src)[static_cast<std::size_t>(lin)];
        // Convert spectral phase [−π,π] to sine-table position [0, SinSiz):
        //   SinBuf[i] = sin(2π·i/SinSiz)  →  phase θ maps to i = θ·SinSiz/(2π)
        double phs = theta / two_pi * SinSiz_d;
        phs = std::fmod(phs, SinSiz_d);
        if (phs < 0.0) phs += SinSiz_d;
        // SAIL: OscPhs is INTEGER ARRAY — truncate (Parshl-source.txt L1046, L1667).
        // [§68] T4 S3: if t4_phase_coherent, skip integer truncation for full phase precision.
        if (t4_phase_coherent) {
          synth_state.OscPhs[osc_sz] = phs;  // [§68] S3: full double precision
        } else {
          synth_state.OscPhs[osc_sz] = static_cast<double>(static_cast<int>(phs));  // SAIL faithful
        }
      }
    }
    // --------------------------------------------------------------------------

    // [§68] T4 S2: soft spectral phase correction for continuing oscillators.
    // Applied before synthesis (after birth init which handled new oscillators).
    // For each continuing osc (PrvOscAmp > 0 and currently active):
    //   phi_synth = OscPhs / SinSiz * 2π   (synthesiser's accumulated phase, radians)
    //   phi_obs   = LinPhase[partial]       (spectral phase from STFT, radians)
    //   delta     = wrap(phi_obs - phi_synth) ∈ [−π, π]
    //   OscPhs   += alpha * delta / (2π) * SinSiz   (corrected, no int truncation)
    if (t4_phase_coherent && t4_alpha > 0.0 && synth_enabled && !reverse_analysis) {
      constexpr double two_pi_t4 = 2.0 * std::numbers::pi;
      const double SinSiz_t4 = static_cast<double>(parshl::ParshlSynthState::SinSiz);
      const auto& st_t4 = tracker.state();
      const std::int64_t lim_t4 = std::min(
          static_cast<std::int64_t>(st_t4.Noscs),
          static_cast<std::int64_t>(synth_state.OscPhs.size()) - 1);
      std::int64_t t4_corr_this_frame = 0;
      for (std::int64_t osc = 1; osc <= lim_t4; ++osc) {
        const std::size_t oi = static_cast<std::size_t>(osc);
        if (oi >= st_t4.LinOfOsc.size() || st_t4.LinOfOsc[oi] <= 0) continue;
        // Only continuing oscillators (birth = handled by S3 above)
        if (oi >= st_t4.PrvOscAmp.size() || st_t4.PrvOscAmp[oi] == 0.0) continue;
        const std::int64_t pk = st_t4.LinOfOsc[oi];
        if (pk < 1 || static_cast<std::size_t>(pk) >= LinPhase.size()) continue;
        // Convert current OscPhs (table units) to radians
        const double phi_synth = synth_state.OscPhs[oi] / SinSiz_t4 * two_pi_t4;
        const double phi_obs   = LinPhase[static_cast<std::size_t>(pk)];
        // Wrapped angular error: delta ∈ [−π, π]
        double delta = std::fmod(phi_obs - phi_synth + std::numbers::pi, two_pi_t4);
        if (delta < 0.0) delta += two_pi_t4;
        delta -= std::numbers::pi;
        // Apply correction: phi_new = phi_synth + alpha * delta
        const double phi_new = phi_synth + t4_alpha * delta;
        double phs_new = phi_new / two_pi_t4 * SinSiz_t4;
        phs_new = std::fmod(phs_new, SinSiz_t4);
        if (phs_new < 0.0) phs_new += SinSiz_t4;
        synth_state.OscPhs[oi] = phs_new;  // no int truncation
        ++t4_corr_this_frame;
      }
      t4_correction_count_accum += t4_corr_this_frame;
    }
    // --------------------------------------------------------------------------

    // [§65] phase_error: McAulay-Quatieri spectral phase continuity error.
    // For each matched oscillator with phase history: |phi_t - (phi_{t-1} + omega*dt)|,
    // wrapped to [0, pi].  phi from atan2(Im,Re) at each partial's fractional bin.
    // Measures how well the linear-FM phase predictor tracks the actual spectral phase.
    // Works with or without synthesis enabled; requires f >= 1 (previous frame exists).
    if (!reverse_analysis) {
      const auto& st_ph = tracker.state();
      const std::int64_t lim_ph = std::min(maxoscs_val,
          static_cast<std::int64_t>(st_ph.LinOfOsc.size()) - 1);
      if (diag_prv_phase.size() < static_cast<std::size_t>(lim_ph + 1))
        diag_prv_phase.assign(static_cast<std::size_t>(lim_ph + 1), 0.0);
      constexpr double two_pi = 2.0 * std::numbers::pi;
      for (std::int64_t osc = 1; osc <= lim_ph; ++osc) {
        const std::size_t oi = static_cast<std::size_t>(osc);
        if (oi >= st_ph.LinOfOsc.size() || st_ph.LinOfOsc[oi] <= 0) continue;
        const std::int64_t p = st_ph.LinOfOsc[oi];
        if (p < 1 || static_cast<std::size_t>(p) >= LinPhase.size()) continue;
        const double phi_t = LinPhase[static_cast<std::size_t>(p)];
        // Measure phase prediction error only for oscillators with phase history
        if (oi < st_ph.PrvOscAmp.size() && st_ph.PrvOscAmp[oi] != 0.0) {
          const double omega   = two_pi * st_ph.OscFrq[oi];
          const double phi_exp = diag_prv_phase[oi]
                                 + omega * static_cast<double>(hop)
                                         / static_cast<double>(fs);
          double pe = std::fmod(std::fabs(phi_t - phi_exp), two_pi);
          if (pe > std::numbers::pi) pe = two_pi - pe;
          diag_phase_err_sum += pe;
          diag_phase_err_max  = std::max(diag_phase_err_max, pe);
          ++diag_phase_err_n;
        }
        // Always update stored phase for next frame
        diag_prv_phase[oi] = phi_t;
      }
    }

    if (flt_enabled) {
      for (std::int64_t i = 1; i <= Nspec; ++i) {
        const double xr = stft.lastS[2*i - 1];
        const double xi = stft.lastS[2*i];
        const double hr = flt_Hr[static_cast<std::size_t>(i - 1)];
        const double hi = flt_Hi[static_cast<std::size_t>(i - 1)];
        ifft.in[static_cast<std::size_t>(i - 1)][0] = xr*hr - xi*hi;
        ifft.in[static_cast<std::size_t>(i - 1)][1] = xr*hi + xi*hr;
      }
      ifft.exec();
      for (int i = 1; i <= synth_nfft; ++i) {
        OutBuf[static_cast<std::size_t>(Bp + i - 1)] += static_cast<float>(ifft.out[static_cast<std::size_t>(i - 1)]);
      }
    }

    if (synth_enabled) {
      // D4: use pre-computed reversed-then-flipped frame parameters when active.
      // Otherwise fall back to the live tracker state (existing behaviour).
      const bool use_rev = reverse_analysis
                           && static_cast<std::size_t>(f) < rev_frames.size();
      const auto& sst = tracker.state(); // live (always available for logging)

      // Pointers to the amp/frq arrays used for synthesis.
      const std::vector<double>* pa_OscAmp    = &sst.OscAmp;
      const std::vector<double>* pa_OscFrq    = &sst.OscFrq;
      const std::vector<double>* pa_PrvOscAmp = &sst.PrvOscAmp;
      const std::vector<double>* pa_PrvOscFrq = &sst.PrvOscFrq;
      std::int64_t synth_Noscs = sst.Noscs;

      if (use_rev) {
        const auto& rf = rev_frames[static_cast<std::size_t>(f)];
        pa_OscAmp    = &rf.OscAmp;
        pa_OscFrq    = &rf.OscFrq;
        pa_PrvOscAmp = &rf.PrvOscAmp;
        pa_PrvOscFrq = &rf.PrvOscFrq;
        synth_Noscs  = rf.Noscs;
      }

      // SAIL: Synthesize(OscAmp, OscFrq, PrvOscAmp, PrvOscFrq, ...) — additive OLA synthesis
      parshl::parshl_synthesize_additive(
        synth_state,
        static_cast<int>(hop),
        Bp,
        fs,
        static_cast<int>(f + 1),
        static_cast<int>(use_rev ? synth_Noscs : tp.MaxOscs),
        *pa_PrvOscAmp,
        *pa_PrvOscFrq,
        *pa_OscAmp,
        *pa_OscFrq,
        synth_skip,
        OutBuf,
        oscphs_double || t4_phase_coherent  // [V3 E4] / [§68] T4 implies full phase precision
      );

      if (synth_log) {
        double e = 0.0;
        double peak = 0.0;
        const std::size_t beg = static_cast<std::size_t>(std::max(0, Bp));
        const std::size_t end = static_cast<std::size_t>(Bp + hop);
        const std::size_t lim = std::min(end, OutBuf.size());
        for (std::size_t i = beg; i < lim; ++i) {
          const double v = static_cast<double>(OutBuf[i]);
          e += v * v;
          peak = std::max(peak, std::abs(v));
        }
        const double rms = (lim > beg) ? std::sqrt(e / static_cast<double>(lim - beg)) : 0.0;
        std::cout << "[SYNTH] frame1=" << (f + 1)
                  << " Bp=" << Bp
                  << " Nhop=" << hop
                  << " Noscs=" << synth_Noscs
                  << " rms=" << rms
                  << " peak=" << peak
                  << (use_rev ? " [D4-rev]" : "")
                  << "\n";
      }

      // Accumulate synthesis stats (peak_abs, overshoot, RMS) from the hot OLA region.
      {
        const std::size_t beg = static_cast<std::size_t>(std::max(0, Bp));
        const std::size_t end = static_cast<std::size_t>(Bp + hop);
        const std::size_t lim = std::min(end, OutBuf.size());
        for (std::size_t i = beg; i < lim; ++i) {
          const double v = static_cast<double>(OutBuf[i]);
          const double av = std::abs(v);
          if (av > taccum.peak_abs) taccum.peak_abs = av;
          if (av > 1.0) ++taccum.overshoot_count;
          taccum.sum_sq += v * v;
          ++taccum.n_samples;
        }
      }

    }

    if (out_enabled) {
      // SAIL: Bp += Nhop; IF Bp >= Nfft THEN flush+decimate OLA buffer, Bp -= Nfft
      // Matches the SAIL "IF Bp >= Nfft THEN" branch in the main loop OLA section.
      Bp += static_cast<int>(hop);
      if (Bp >= synth_nfft) {
        const int I = decim_first_call ? 1 : 0;
        const int Nout = decimate_sail(OutBuf, synth_nfft, Ndec, I, decim_P);
        decim_first_call = false;
        out_buf.insert(out_buf.end(), OutBuf.begin(), OutBuf.begin() + Nout);
        std::copy(OutBuf.begin() + synth_nfft, OutBuf.begin() + 2 * synth_nfft, OutBuf.begin());
        std::fill(OutBuf.begin() + synth_nfft, OutBuf.begin() + 2 * synth_nfft, 0.0f);
        Bp -= synth_nfft;
      }
      out_written = out_buf.size() + static_cast<std::size_t>(std::max(0, Bp));
    }

    // ---- FindPeaks (apples-to-apples comparison) ----
    // IMPORTANT: FindPeaks clobbers X, so we copy the dB spectrum.
    sail::Array1D<double> Xtmp = Xdb;

    const std::int64_t BinInt = BinInt_from_MinSepHz(MinSepHz, (double)fs, Nfft);

    const int npeaks = parshl::FindPeaks(
      Xtmp,
      Peaks,
      Locs,
      thresh_db,
      Hyst,
      /*MaxPeaks=*/tp.MaxLins,
      /*MinWidth=*/MinWid,
      /*BlankWidth=*/BinInt,
      /*I1=*/I1,
      /*I2=*/I2
    );

    std::cout << "Frame " << f << " pos=" << pos << " peaks=" << npeaks << "\n";
    for (int i = 1; i <= npeaks && i <= 10; ++i) {
      std::cout << "  " << i << ": dB=" << Peaks[i] << " bin=" << Locs[i] << "\n";
    }

    // ---- InterpComplexAtPeak (debug only) ----
    // Paper: Serra & Smith, ICMC 1990 §5 Step 4 — NO SAIL EQUIVALENT.
    // NOT connected to the main analysis pipeline (Qinterp on magnitude only).
    // See docs/COMPLEX_PEAK_AUDIT.md for full analysis.
    const int toShow = std::min(npeaks, 5);
    for (int i = 1; i <= toShow; ++i) {
      const double loc = Locs[i];

      // k = nearest bin (1-based)
      std::int64_t k = static_cast<std::int64_t>(std::llround(loc));
      double p = loc - static_cast<double>(k);

      if (k < 1) k = 1;
      if (k > Nspec) k = Nspec;

      double re = 0.0, im = 0.0;
      parshl::InterpComplexAtPeak(stft.lastS, (long)Nspec, (long)k, p, re, im);

      const double mag = std::sqrt(re*re + im*im);
      const double ph  = std::atan2(im, re);

      std::cout << "    complex@" << i
                << ": k=" << k << " p=" << p
                << " re=" << re << " im=" << im
                << " mag=" << mag << " phase=" << ph << "\n";
    }

    // ---- Tracks ----
    const auto& st = tracker.state();
    std::cout << "  Tracks: Noscs=" << st.Noscs;
    if (reverse_analysis && static_cast<std::size_t>(f) < rev_frames.size()) {
      std::cout << " [D4-rev Noscs=" << rev_frames[static_cast<std::size_t>(f)].Noscs << "]";
    }
    std::cout << "\n";
    for (std::int64_t osc = 1; osc <= st.Noscs && osc <= 10; ++osc) {
      std::cout << "    osc " << osc
                << " lin=" << st.LinOfOsc[osc]
                << " amp=" << st.OscAmp[osc]
                << " frqHz=" << st.OscFrq[osc] << "\n";
    }

    // ---- R1 OutPartials write (§28) ----
    // Writes OscAmp and OscFrq for every oscillator in the bank to the trajectory WAV files.
    // For D4 (reverse_analysis): use the pre-computed rev_frames parameters (same source
    // as synthesis), so the trajectories reflect the time-corrected forward-order amplitudes.
    // For modes A/B: use the live tracker state directly.
    // This block is a pure side-output — zero effect on tracker state, synthesis, or stats.
    if (partials_writer.is_open()) {
      const bool use_rev_pw = reverse_analysis
                              && static_cast<std::size_t>(f) < rev_frames.size();
      const std::vector<double>& pw_amp = use_rev_pw
          ? rev_frames[static_cast<std::size_t>(f)].OscAmp
          : tracker.state().OscAmp;
      const std::vector<double>& pw_frq = use_rev_pw
          ? rev_frames[static_cast<std::size_t>(f)].OscFrq
          : tracker.state().OscFrq;
      partials_writer.write_frame(pw_amp, pw_frq, maxoscs_val);
    }

    pos += hop;
    if (pos >= (std::int64_t)x.size()) break;
  }

  if (out_enabled) {
    // SAIL final flush uses Decimate(OutBuf, Bp+Nh, ...).
    const int Nh = flt_enabled ? flt_Nh : 0;
    const int tail_len = Bp + Nh;
    if (tail_len > 0) {
      const int I = decim_first_call ? 1 : 0;
      const int Nout = decimate_sail(OutBuf, tail_len, Ndec, I, decim_P);
      decim_first_call = false;
      out_buf.insert(out_buf.end(), OutBuf.begin(), OutBuf.begin() + Nout);
    }
    out_written = out_buf.size();
    std::vector<double> out_d(out_written, 0.0);
    for (std::size_t i = 0; i < out_written; ++i) out_d[i] = static_cast<double>(out_buf[i]);
    // [FIX] Normalise from SAIL units to [-1.0, 1.0] before writing synthesis WAV.
    // The ×32768 applied on read (io_sndfile.cpp) must be reversed so that the
    // float32 WAV is playable in any audio player.
    for (std::size_t i = 0; i < out_written; ++i) out_d[i] /= 32768.0;
    // [N1] Optional post-synthesis peak-level normalization (--normalize-peak=<dBFS>).
    // Applied after the 32768 rescale; SDR is computed from un-normalized data above.
    if (normalize_peak_enabled) {
      double cur_peak = 0.0;
      for (std::size_t si = 0; si < out_written; ++si)
        cur_peak = std::max(cur_peak, std::abs(out_d[si]));
      if (cur_peak > 0.0) {
        const double target = std::pow(10.0, normalize_peak_dbfs / 20.0);
        const double gain   = target / cur_peak;
        for (std::size_t si = 0; si < out_written; ++si) out_d[si] *= gain;
        std::cout << "[NORM] peak-normalize: cur_peak=" << cur_peak
                  << "  target=" << target
                  << "  gain=" << 20.0 * std::log10(gain) << " dB\n";
      }
    }
    if (!parshl::WriteMonoWavFloat32(synth_out_wav, out_d, fs)) {
      std::cerr << "failed writing synthesized wav: " << synth_out_wav << "\n";
      return 4;
    }
    std::cout << "[SYNTH] wrote " << out_written << " samples to " << synth_out_wav << "\n";
  }

  // ---- Close still-alive trajectories at end of run (modes A/B only) -----------
  // Mode C trajectory accounting was done in the pre-loop compute_revframes_lifecycle().
  if (!reverse_analysis) {
    const auto& st_final = tracker.state();
    const std::int64_t lim_f = std::min(maxoscs_val,
                                        static_cast<std::int64_t>(st_final.LinOfOsc.size()) - 1);
    for (std::int64_t osc = 1; osc <= lim_f; ++osc) {
      const std::size_t oi = static_cast<std::size_t>(osc);
      if (osc_birth[oi] >= 0 && oi < st_final.LinOfOsc.size()) {
        // trajectory alive at last frame — close with remaining duration
        taccum.lifetimes.push_back(
            static_cast<std::int64_t>(taccum.n_frames) - osc_birth[oi]);
      }
    }
    // Sync add/stop/recycle/dfmax from tracker.stats() (authoritative source).
    taccum.add_count     = tracker.stats().add_count;
    taccum.stop_count    = tracker.stats().stop_count;
    taccum.recycle_count = tracker.stats().recycle_count;
  }

  // ---- Compute summary metrics -------------------------------------------------
  const std::int64_t track_fragments =
      taccum.add_count + taccum.recycle_count;  // each start = one fragment

  double avg_lifetime = 0.0;
  double median_lifetime = 0.0;
  if (!taccum.lifetimes.empty()) {
    double sum_lt = 0.0;
    for (auto lt : taccum.lifetimes) sum_lt += static_cast<double>(lt);
    avg_lifetime = sum_lt / static_cast<double>(taccum.lifetimes.size());

    std::vector<std::int64_t> sorted_lt = taccum.lifetimes;
    std::sort(sorted_lt.begin(), sorted_lt.end());
    const std::size_t mid = sorted_lt.size() / 2;
    median_lifetime = (sorted_lt.size() % 2 == 0)
        ? (static_cast<double>(sorted_lt[mid - 1]) + static_cast<double>(sorted_lt[mid])) / 2.0
        : static_cast<double>(sorted_lt[mid]);
  }

  const double avg_partials =
      (taccum.n_frames > 0) ? (static_cast<double>(taccum.total_nlins) / static_cast<double>(taccum.n_frames)) : 0.0;
  const double avg_noscs =
      (taccum.n_frames > 0) ? (static_cast<double>(taccum.total_noscs) / static_cast<double>(taccum.n_frames)) : 0.0;
  const double stop_add_ratio =
      (track_fragments > 0) ? (static_cast<double>(taccum.stop_count) / static_cast<double>(track_fragments)) : 0.0;
  const double rms_val =
      (taccum.n_samples > 0) ? std::sqrt(taccum.sum_sq / static_cast<double>(taccum.n_samples)) : 0.0;

  const std::int64_t dfmax_ex = reverse_analysis
      ? rev_dfmax_exceeded                  // from tracker_rev (reversed pass)
      : tracker.stats().dfmax_exceeded;     // from live forward tracker (modes A/B)

  // [§32B] DBspec refinement: compute delta distribution (mean, median, max).
  double dbspec_mean_delta   = 0.0;
  double dbspec_median_delta = 0.0;
  double dbspec_max_delta    = 0.0;
  if (!dbspec_deltas.empty()) {
    for (double d : dbspec_deltas) {
      dbspec_mean_delta += d;
      dbspec_max_delta   = std::max(dbspec_max_delta, d);
    }
    dbspec_mean_delta /= static_cast<double>(dbspec_deltas.size());
    std::vector<double> sorted_dd = dbspec_deltas;
    std::sort(sorted_dd.begin(), sorted_dd.end());
    const std::size_t mid_dd = sorted_dd.size() / 2;
    dbspec_median_delta = (sorted_dd.size() % 2 == 0)
        ? (sorted_dd[mid_dd - 1] + sorted_dd[mid_dd]) / 2.0
        : sorted_dd[mid_dd];
  }

  // [V3 E3] DBspec amplitude delta distribution (mean, median, max in linear amplitude units).
  double dbspec_amp_mean_delta   = 0.0;
  double dbspec_amp_median_delta = 0.0;
  double dbspec_amp_max_delta    = 0.0;
  if (!dbspec_amp_deltas.empty()) {
    for (double d : dbspec_amp_deltas) {
      dbspec_amp_mean_delta += d;
      dbspec_amp_max_delta   = std::max(dbspec_amp_max_delta, d);
    }
    dbspec_amp_mean_delta /= static_cast<double>(dbspec_amp_deltas.size());
    std::vector<double> sorted_ad = dbspec_amp_deltas;
    std::sort(sorted_ad.begin(), sorted_ad.end());
    const std::size_t mid_ad = sorted_ad.size() / 2;
    dbspec_amp_median_delta = (sorted_ad.size() % 2 == 0)
        ? (sorted_ad[mid_ad - 1] + sorted_ad[mid_ad]) / 2.0
        : sorted_ad[mid_ad];
  }

  // ---- Print STATS_SUMMARY (parseable by extract_metrics.py) ------------------
  std::cout << std::fixed << std::setprecision(6)
            << "[STATS_SUMMARY]"
            << " n_frames="        << taccum.n_frames
            << " avg_partials="    << avg_partials
            << " max_partials="    << taccum.max_nlins
            << " Noscs_avg="       << avg_noscs
            << " Noscs_max="       << taccum.max_noscs
            << " add_count="       << taccum.add_count
            << " stop_count="      << taccum.stop_count
            << " recycle_count="   << taccum.recycle_count
            << " track_fragments=" << track_fragments
            << " avg_lifetime="    << avg_lifetime
            << " median_lifetime=" << median_lifetime
            << " dfmax_exceeded="  << dfmax_ex
            << " damax_exceeded_observe=" << tracker.stats().damax_exceeded_observe
            << " damax_gate_stops="       << tracker.stats().damax_gate_stops
            << " osc_recycles="    << taccum.recycle_count
            << " peak_abs="        << taccum.peak_abs
            << " overshoot_count=" << taccum.overshoot_count
            << " rms="             << rms_val
            << " reverse_analysis=" << (reverse_analysis ? 1 : 0)
            << " legacy="          << (legacy_mode ? 1 : 0)
            << " damax_observe="        << tp.DAmax_observe
            << " damax_gate="           << tp.DAmax_gate
            << " dbspec_refine="        << (dbspec_refine ? 1 : 0)
            << " dbspec_n_refined="     << dbspec_deltas.size()
            << " dbspec_mean_delta_hz=" << dbspec_mean_delta
            << " dbspec_median_delta_hz=" << dbspec_median_delta
            << " dbspec_max_delta_hz="  << dbspec_max_delta
            << " closest_osc_free="     << (closest_osc_free ? 1 : 0)
            << " cof_reuse_count="      << tracker.stats().cof_reuse_count
            << " sigscl_analytic="      << (sigscl_analytic ? 1 : 0)
            << " sigscl_value="         << std::setprecision(10) << SigScl
            << " sigscl_sail_value="    << std::setprecision(10) << sail_sigscl
            << " dbspec_refine_amp="        << (dbspec_refine_amp ? 1 : 0)
            << " dbspec_amp_n_refined="     << dbspec_amp_deltas.size()
            << " dbspec_amp_mean_delta="    << std::setprecision(10) << dbspec_amp_mean_delta
            << " dbspec_amp_median_delta="  << std::setprecision(10) << dbspec_amp_median_delta
            << " dbspec_amp_max_delta="     << std::setprecision(10) << dbspec_amp_max_delta
            << " oscphs_double="             << (oscphs_double ? 1 : 0)
            << " dbspec_match="             << (dbspec_match ? 1 : 0)
            << " dbspec_match_count="       << dbspec_match_count
            << " grs_groups="                << grs_groups
            << " grs_group_collision_count=" << grs_group_collision_count
            << " t1_peak_claim="             << (t1_peak_claim ? 1 : 0)
            << " t1_created="               << t1_stats_accum.t1_created
            << " t1_assigned="              << t1_stats_accum.t1_assigned
            << " t1_terminated="            << t1_stats_accum.t1_terminated
            << " t1_squelch_hold="          << t1_squelch_hold
            << " t1_hold_survivals="        << t1_stats_accum.t1_hold_survivals
            << " t1_hold_terminations="     << t1_stats_accum.t1_hold_terminations
            << " t1_amp_hysteresis="        << t1_amp_hysteresis  // [§64]
            << " t2_global_assign="         << (t2_global_assign ? 1 : 0)
            << " t2_assignments="           << t2_stats_accum.assignments
            << " t2_creations="             << t2_stats_accum.creations
            << " t2_terminations="          << t2_stats_accum.terminations
            << " t2_trajectory_cost="       << (t2_trajectory_cost ? 1 : 0)  // [§60]
            << " t2_identity_switches="     << t2_stats_accum.identity_switches  // [§60]
            << " t2_rank_cost="             << (t2_rank_cost ? 1 : 0)  // [§61]
            << " t2_wrank="                 << t2_wrank                 // [§61]
            // [§66B] Phase-aware Hungarian cost
            << " t2_phase_cost="            << (t2_phase_cost ? 1 : 0)
            << " t2_wphi="                  << t2_wphi
            << " t2_phase_eligible="        << t2_stats_accum.phase_eligible
            << " t2_mean_phase_cost_rad="   << (t2_stats_accum.phase_eligible > 0
                ? t2_stats_accum.phase_cost_sum / static_cast<double>(t2_stats_accum.phase_eligible)
                : 0.0)
            // [§68] T4 phase-coherent synthesis
            << " t4_phase_coherent="          << (t4_phase_coherent ? 1 : 0)
            << " t4_alpha="                   << t4_alpha
            << " t4_corrections="             << t4_correction_count_accum
            // [§65] Analysis error diagnostics
            << " diag_mean_freq_err_hz="  << (diag_freq_err_n  > 0
                ? diag_freq_err_sum  / static_cast<double>(diag_freq_err_n)  : 0.0)
            << " diag_max_freq_err_hz="   << diag_freq_err_max
            << " diag_mean_amp_err="      << (diag_amp_err_n   > 0
                ? diag_amp_err_sum   / static_cast<double>(diag_amp_err_n)   : 0.0)
            << " diag_mean_phase_err_rad=" << (diag_phase_err_n > 0
                ? diag_phase_err_sum / static_cast<double>(diag_phase_err_n) : 0.0)
            << " diag_max_phase_err_rad=" << diag_phase_err_max
            << "\n";

  return 0;
}
