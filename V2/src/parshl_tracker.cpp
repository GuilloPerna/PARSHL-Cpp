// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "parshl_tracker.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <stdexcept>

namespace parshl {

namespace {

static inline std::int64_t env_i64(const char* name, std::int64_t defv) {
  const char* s = std::getenv(name);
  if (!s || !*s) return defv;
  char* end = nullptr;
  const long long v = std::strtoll(s, &end, 10);
  if (end == s) return defv;
  return static_cast<std::int64_t>(v);
}

static inline bool trk_log_enabled() {
  static const bool on = (env_i64("PARSHL_TRK_LOG", 0) != 0);
  return on;
}

static inline bool trk_checks_enabled() {
  static const bool on = (env_i64("PARSHL_TRK_CHECKS", 0) != 0);
  return on;
}

static inline bool trk_snapshot_enabled() {
  static const bool on = (env_i64("PARSHL_TRK_SNAPSHOT", 0) != 0);
  return on;
}

static inline bool trk_trace_writes_enabled() {
  static const bool on = (env_i64("PARSHL_TRK_TRACE_WRITES", 0) != 0);
  return on;
}

static inline bool trk_dfmax_stats_enabled() {
  static const bool on = (env_i64("PARSHL_TRK_DFMAX_STATS", 0) != 0);
  return on;
}

static inline bool trk_shadow_enabled() {
  static const bool on = (env_i64("PARSHL_TRK_SHADOW", 0) != 0);
  return on;
}

static inline bool trk_verbose_enabled() {
  static const bool on = (env_i64("PARSHL_TRK_VERBOSE", 0) != 0);
  return on;
}

static inline std::int64_t trk_frame0() {
  static const std::int64_t v = env_i64("PARSHL_TRK_FRAME0", 0);
  return v;
}

static inline std::int64_t trk_frame1_limit() {
  static const std::int64_t v = env_i64("PARSHL_TRK_FRAME1", -1);
  return v;
}

static inline std::int64_t trk_osc_max() {
  static const std::int64_t v = env_i64("PARSHL_TRK_OSC_MAX", 0);
  return v;
}

static inline std::int64_t trk_lin_max() {
  static const std::int64_t v = env_i64("PARSHL_TRK_LIN_MAX", 0);
  return v;
}

static inline std::int64_t trk_case_max() {
  static const std::int64_t v = env_i64("PARSHL_TRK_CASE_MAX", 12);
  return (v > 0) ? v : 12;
}

static inline std::int64_t trk_topk() {
  static const std::int64_t v = env_i64("PARSHL_TRK_TOPK", 10);
  return (v > 0) ? v : 10;
}

struct WriteMeta {
  std::int64_t frame1 = -1;
  const char* func = "unset";
  const char* reason = "unset";
  std::int64_t oldv = 0;
  std::int64_t newv = 0;
  std::int64_t seq = 0;
};

struct RestoreEvent {
  std::int64_t lin = 0;
  std::int64_t savedOsc = 0;
  std::int64_t curOsc = 0;
  std::int64_t otherOsc = 0;
};

struct RatioRec {
  double ratio = 0.0;
  std::int64_t osc = 0;
  std::int64_t bestLin = 0;
  double PF = 0.0;
  double BestFrq = 0.0;
  double MinDist = 0.0;
  double DFmax = 0.0;
  const char* decision = "unknown";
};

thread_local std::int64_t g_frame1 = -1;
thread_local const char* g_stop_reason = "unset";
thread_local int g_dfmax_kase = 0;
thread_local int g_rec_depth = 0;

thread_local std::int64_t g_dfmax_exceeded_frame = 0;
thread_local std::int64_t g_no_best_line_frame = 0;
thread_local std::int64_t g_collision_count_frame = 0;
thread_local std::int64_t g_recursion_depth_max_frame = 0;
thread_local std::int64_t g_invariants_fail_frame = 0;
thread_local std::int64_t g_stop_count_frame = 0;
thread_local std::int64_t g_add_count_frame = 0;
thread_local std::int64_t g_i2_fail_count_frame = 0;
thread_local std::int64_t g_i11_fail_count_frame = 0;
thread_local std::int64_t g_collision_restore_saved_frame = 0;

thread_local std::int64_t g_dfmax_exceeded_acc = 0;
thread_local std::int64_t g_no_best_line_acc = 0;
thread_local std::int64_t g_collision_count_acc = 0;
thread_local std::int64_t g_recursion_depth_max_acc = 0;
thread_local std::int64_t g_invariants_fail_acc = 0;
thread_local std::int64_t g_stop_count_acc = 0;
thread_local std::int64_t g_add_count_acc = 0;
thread_local std::int64_t g_i2_fail_count_acc = 0;
thread_local std::int64_t g_i11_fail_count_acc = 0;
thread_local std::int64_t g_collision_restore_saved_acc = 0;

thread_local std::int64_t g_shadow_unclaim_delta_frame = 0;
thread_local std::int64_t g_shadow_unclaim_delta_acc = 0;
thread_local std::int64_t g_shadow_instantrise_delta_frame = 0;
thread_local std::int64_t g_shadow_instantrise_delta_acc = 0;

thread_local std::array<std::int64_t, 4> g_dfmax_ratio_hist_frame{{0, 0, 0, 0}};
thread_local std::array<std::int64_t, 4> g_dfmax_ratio_hist_acc{{0, 0, 0, 0}};
thread_local std::vector<RatioRec> g_dfmax_top_frame;

thread_local std::vector<WriteMeta> g_lastwrite_linofosc;
thread_local std::vector<WriteMeta> g_lastwrite_oscoflin;
thread_local std::int64_t g_write_seq = 0;

thread_local std::vector<RestoreEvent> g_restore_events_frame;

thread_local std::vector<int> g_osc_action; // 0=none, 1=take, 2=stop, 3=add

static inline bool trk_frame_in_window() {
  if (!(trk_log_enabled() || trk_checks_enabled() || trk_snapshot_enabled() ||
        trk_trace_writes_enabled() || trk_dfmax_stats_enabled() || trk_shadow_enabled())) {
    return false;
  }
  if (g_frame1 < 1) return true;
  const std::int64_t frame0 = g_frame1 - 1;
  const std::int64_t f0 = trk_frame0();
  const std::int64_t f1 = trk_frame1_limit();
  if (f1 < 0) return (frame0 >= f0);
  return (frame0 >= f0 && frame0 <= f1);
}

static inline bool trk_log_osc(std::int64_t osc) {
  if (!trk_frame_in_window()) return false;
  if (osc < 1) return false;
  const std::int64_t lim = trk_osc_max();
  return (lim <= 0) ? true : (osc <= lim);
}

static inline void trk_log(const std::string& msg) {
  if (!trk_frame_in_window()) return;
  std::cerr << "[TRK] " << msg << "\n";
}

static inline std::uint64_t hash_mix(std::uint64_t h, std::uint64_t x) {
  h ^= x;
  h *= 1099511628211ULL;
  return h;
}

static inline std::uint64_t hash_i64_vec(const std::vector<std::int64_t>& v) {
  std::uint64_t h = 1469598103934665603ULL;
  for (std::size_t i = 1; i < v.size(); ++i) {
    h = hash_mix(h, static_cast<std::uint64_t>(v[i]));
  }
  return h;
}

static inline std::uint64_t hash_f64_vec(const std::vector<double>& v) {
  std::uint64_t h = 1469598103934665603ULL;
  for (std::size_t i = 1; i < v.size(); ++i) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &v[i], sizeof(bits));
    h = hash_mix(h, bits);
  }
  return h;
}

static inline void ensure_write_meta_sizes(std::int64_t maxOscs, std::int64_t maxLins) {
  const std::size_t osz = static_cast<std::size_t>(maxOscs + 1);
  const std::size_t lsz = static_cast<std::size_t>(maxLins + 1);
  if (g_lastwrite_linofosc.size() != osz) g_lastwrite_linofosc.assign(osz, WriteMeta{});
  if (g_lastwrite_oscoflin.size() != lsz) g_lastwrite_oscoflin.assign(lsz, WriteMeta{});
}

static inline void mark_linofosc_write(std::int64_t osc, std::int64_t oldv, std::int64_t newv,
                                       const char* func, const char* reason) {
  if (!trk_trace_writes_enabled()) return;
  if (osc <= 0 || osc >= static_cast<std::int64_t>(g_lastwrite_linofosc.size())) return;
  WriteMeta m;
  m.frame1 = g_frame1;
  m.func = func;
  m.reason = reason;
  m.oldv = oldv;
  m.newv = newv;
  m.seq = ++g_write_seq;
  g_lastwrite_linofosc[static_cast<std::size_t>(osc)] = m;
}

static inline void mark_oscoflin_write(std::int64_t lin, std::int64_t oldv, std::int64_t newv,
                                       const char* func, const char* reason) {
  if (!trk_trace_writes_enabled()) return;
  if (lin <= 0 || lin >= static_cast<std::int64_t>(g_lastwrite_oscoflin.size())) return;
  WriteMeta m;
  m.frame1 = g_frame1;
  m.func = func;
  m.reason = reason;
  m.oldv = oldv;
  m.newv = newv;
  m.seq = ++g_write_seq;
  g_lastwrite_oscoflin[static_cast<std::size_t>(lin)] = m;
}

static inline std::string fmt_write_meta(const WriteMeta& m) {
  std::ostringstream oss;
  oss << "{frame1=" << m.frame1
      << " seq=" << m.seq
      << " func=" << m.func
      << " reason=" << m.reason
      << " old=" << m.oldv
      << " new=" << m.newv
      << "}";
  return oss.str();
}

static inline int ratio_bin(double r) {
  if (r <= 1.0) return 0;
  if (r <= 2.0) return 1;
  if (r <= 5.0) return 2;
  return 3;
}

static inline void record_ratio(const RatioRec& rec) {
  if (!trk_dfmax_stats_enabled()) return;
  const int b = ratio_bin(rec.ratio);
  ++g_dfmax_ratio_hist_frame[static_cast<std::size_t>(b)];
  ++g_dfmax_ratio_hist_acc[static_cast<std::size_t>(b)];

  g_dfmax_top_frame.push_back(rec);
  std::sort(g_dfmax_top_frame.begin(), g_dfmax_top_frame.end(),
            [](const RatioRec& a, const RatioRec& b) { return a.ratio > b.ratio; });
  const std::size_t k = static_cast<std::size_t>(trk_topk());
  if (g_dfmax_top_frame.size() > k) g_dfmax_top_frame.resize(k);
}

static inline std::string dump_i64_1based(const std::vector<std::int64_t>& v, std::int64_t k) {
  std::ostringstream oss;
  oss << "[";
  for (std::int64_t i = 1; i <= k; ++i) {
    if (i > 1) oss << " ";
    if (i < static_cast<std::int64_t>(v.size())) oss << i << ":" << v[static_cast<std::size_t>(i)];
    else oss << i << ":<oob>";
  }
  oss << "]";
  return oss.str();
}

} // namespace

static inline double absd(double x) { return std::fabs(x); }

// Access helpers that tolerate either 1-based (size >= Nlins+1) or 0-based (size == Nlins)
// input vectors. The header says 1..Nlins, but this avoids errors if the caller passed 0-based.
static inline double get_lin(const std::vector<double>& v, std::int64_t idx_1based, std::int64_t Nlins) {
  if (idx_1based < 1 || idx_1based > Nlins) return 0.0;
  const auto n = static_cast<std::int64_t>(v.size());
  if (n >= Nlins + 1) return v[static_cast<std::size_t>(idx_1based)];
  if (n == Nlins)     return v[static_cast<std::size_t>(idx_1based - 1)];
  // fallback: best effort
  if (!v.empty()) {
    const auto z = std::min<std::size_t>(static_cast<std::size_t>(std::max<std::int64_t>(0, idx_1based - 1)),
                                         v.size() - 1);
    return v[z];
  }
  return 0.0;
}

ParshlTracker::ParshlTracker(const TrackerParams& p) : p_(p) {
  reset();
}

void ParshlTracker::reset() {
  if (p_.MaxOscs <= 0 || p_.MaxLins <= 0) {
    throw std::runtime_error("ParshlTracker: MaxOscs/MaxLins must be > 0");
  }

  st_.Noscs = 0;

  const std::size_t osz = static_cast<std::size_t>(p_.MaxOscs + 1);
  const std::size_t lsz = static_cast<std::size_t>(p_.MaxLins + 1);

  st_.LinOfOsc.assign(osz, 0);
  st_.PrvLinOfOsc.assign(osz, 0);
  st_.OscAmp.assign(osz, 0.0);
  st_.PrvOscAmp.assign(osz, 0.0);
  st_.OscFrq.assign(osz, 0.0);
  st_.PrvOscFrq.assign(osz, 0.0);

  st_.OscOfLin.assign(lsz, 0);

  st_.df_ready = false;
  st_.alpha = 0.0;
  st_.beta = 0.0;

  ts_ = TrackerStats{};
}

// SAIL source: oscillator lifecycle helper
// Parshl-source.txt (within UpdateMap / tracker logic)
//
// Marks oscillator `osc` as Squelched: LinOfOsc[osc] <- -ABS(PrvLinOfOsc[osc])
// The negative sign is the SAIL convention for a decaying (squelched) oscillator.
// A free oscillator (PrvLinOfOsc==0) is returned to free state (LinOfOsc=0).
void ParshlTracker::stop_osc(std::int64_t osc) {
  if (osc < 1 || osc > p_.MaxOscs) return;

  ++g_stop_count_frame;
  ++g_stop_count_acc;
  if (osc < static_cast<std::int64_t>(g_osc_action.size())) {
    g_osc_action[static_cast<std::size_t>(osc)] = 2;
  }

  // PARSHL: LinOfOsc[x] <- -ABS(PrvLinOfOsc[x])
  const std::int64_t prev = (osc < static_cast<std::int64_t>(st_.PrvLinOfOsc.size()))
                              ? st_.PrvLinOfOsc[static_cast<std::size_t>(osc)]
                              : 0;

  // TrackerStats.stop_count: count only TRUE ON→OFF transitions (prev > 0 → cur ≤ 0).
  // Noop stops (prev == 0) and squelch continuations (prev < 0) are excluded so that
  // the count is consistent with the state-transition-based lifecycle counting in the
  // driver (compute_revframes_lifecycle / main-loop lifecycle accumulator).
  if (prev > 0) ++ts_.stop_count;
  if (prev == 0) {
    // If it was free, stopping is effectively a noop -> keep it free
    const std::int64_t oldv = st_.LinOfOsc[static_cast<std::size_t>(osc)];
    st_.LinOfOsc[static_cast<std::size_t>(osc)] = 0;
    mark_linofosc_write(osc, oldv, 0, "stop_osc", g_stop_reason);
  } else {
    const std::int64_t newv = -std::llabs(prev);
    const std::int64_t oldv = st_.LinOfOsc[static_cast<std::size_t>(osc)];
    st_.LinOfOsc[static_cast<std::size_t>(osc)] = newv;
    mark_linofosc_write(osc, oldv, newv, "stop_osc", g_stop_reason);
  }

  if (trk_log_osc(osc)) {
    std::ostringstream oss;
    oss << "stop_osc frame1=" << g_frame1
        << " osc=" << osc
        << " reason=" << g_stop_reason
        << " prvLin=" << prev
        << " newLin=" << st_.LinOfOsc[static_cast<std::size_t>(osc)]
        << " prvFrqHz=" << st_.PrvOscFrq[static_cast<std::size_t>(osc)]
        << " frqHz=" << st_.OscFrq[static_cast<std::size_t>(osc)];
    trk_log(oss.str());
  }
}

// SAIL source: DFmax(f) function — maximum allowed frequency deviation for matching
// Parshl-source.txt (within UpdateMap / GetClosestFrq calling context)
//
// Returns the maximum frequency jump (Hz) allowed when matching oscillator `osc`
// to a new spectral partial.  DFmax is piecewise-linear in frequency:
//   DFmax(f) = alpha * f + beta
// where alpha/beta are computed from the (Fc1, DFmax1) and (Fc2, DFmax2) control points.
// Three cases:
//   kase=1 (general): alpha = (DFmax2/f2 - DFmax1/f1) / (f2-f1)
//   kase=2 (f1=0):    linear pass through origin adjusted by DFmax1
//   kase=3 (constant): alpha=0, beta=DFmax1
double ParshlTracker::dfmax(std::int64_t osc) {
  if (!st_.df_ready) {
    // Reproduce original DFmax init logic (DFmax is linear in frequency).
    // Cases:
    //  - DFmax1 == DFmax2 -> constant
    //  - Fc1 == 0 -> alternate linear form avoiding division by zero
    //  - general case
    const double f1 = std::max(0.0, p_.Fc1);
    const double f2 = std::min(p_.Fc2, p_.Fs * 0.5);

    // Guard degenerate setup
    const double denom = (f2 - f1);
    if (p_.DFmax1 == p_.DFmax2 || denom == 0.0 || f2 <= 0.0) {
      g_dfmax_kase = 3;
      st_.alpha = 0.0;
      st_.beta = p_.DFmax1;
    } else if (f1 == 0.0) {
      // kase = 2
      g_dfmax_kase = 2;
      st_.alpha = (p_.DFmax2 - p_.DFmax1) / denom;
      st_.beta  = (p_.DFmax1 * f2 - p_.DFmax2 * f1) / denom; // with f1=0 -> DFmax1
    } else {
      // kase = 1 (general)
      g_dfmax_kase = 1;
      st_.alpha = ((p_.DFmax2 / f2) - (p_.DFmax1 / f1)) / denom;
      st_.beta  = ((p_.DFmax1 * f2 / f1) - (p_.DFmax2 * f1 / f2)) / denom;
    }

    st_.df_ready = true;
  }

  const double frq = (osc >= 1 && osc <= p_.MaxOscs)
                       ? st_.OscFrq[static_cast<std::size_t>(osc)]
                       : 0.0;
  double out = st_.alpha * frq + st_.beta;

  // V2 D7: adaptive DFmax — expand search window by observed oscillator velocity.
  // velocity = |OscFrq[osc] - PrvOscFrq[osc]| (Hz change over last frame).
  // OutAF guarantees PrvOscFrq == OscFrq when PrvOscAmp == 0, so velocity == 0
  // for newly allocated / recycled oscillators (safe fall-back to SAIL baseline).
  double adaptive_hz = 0.0;
  if (!p_.legacy_mode && p_.dfmax_velocity_scale > 0.0
      && osc >= 1 && osc <= p_.MaxOscs) {
    const std::size_t oi = static_cast<std::size_t>(osc);
    const double velocity = std::fabs(st_.OscFrq[oi] - st_.PrvOscFrq[oi]);
    adaptive_hz = p_.dfmax_velocity_scale * velocity;
    if (adaptive_hz > out) out = adaptive_hz;
  }

  const double fs_half = (p_.Fs > 0.0) ? (0.5 * p_.Fs) : -1.0;
  const double f2_min_fc2_fs2 = (fs_half > 0.0) ? std::min(p_.Fc2, fs_half) : p_.Fc2;
  const double f2_used = std::min(p_.Fc2, p_.Fs * 0.5);

  if (trk_log_osc(osc)) {
    std::ostringstream oss;
    oss << "dfmax frame1=" << g_frame1
      << " osc=" << osc
      << " oscFrqHz=" << frq
      << " alpha=" << st_.alpha
      << " beta=" << st_.beta
      << " kase=" << g_dfmax_kase
      << " dfmaxHz=" << out
      << " adaptiveHz=" << adaptive_hz;
    if (trk_verbose_enabled()) {
      oss << " f2_used_Hz=" << f2_used
        << " Fc2_Hz=" << p_.Fc2
        << " FsOver2_Hz=" << fs_half
        << " minFc2FsOver2_Hz=" << f2_min_fc2_fs2
        << " f2_mode=min_Fc2_FsOver2";
    }
    trk_log(oss.str());

    if (trk_shadow_enabled()) {
      const bool f2_delta = (std::fabs(f2_used - f2_min_fc2_fs2) > 0.0);
      if (f2_delta || trk_verbose_enabled()) {
        std::ostringstream sh;
        sh << "shadow_f2 frame1=" << g_frame1
           << " osc=" << osc
           << " Fc2_Hz=" << p_.Fc2
           << " FsOver2_Hz=" << fs_half
           << " f2_actual_Hz=" << f2_used
           << " f2_shadow_minFc2FsOver2_Hz=" << f2_min_fc2_fs2
           << " deltaHz=" << (f2_used - f2_min_fc2_fs2)
           << " differs=" << (f2_delta ? 1 : 0);
        trk_log(sh.str());
      }
    }
  }

  return out;
}

// SAIL source: PROCEDURE GetClosestFrq (recursive collision resolver)
// Parshl-source.txt (within UpdateMap, step 3)
//
// For oscillator CurOsc, finds the spectral partial (lin) with the closest
// frequency that is within DFmaxVal Hz and is not already claimed.
// If the best candidate is already claimed by another oscillator, recursively
// resolves the collision: the other oscillator is re-assigned or stopped.
// Matches SAIL’s mutual-exclusion assignment algorithm.
//
// D6/DAmax (reconstructive): LinAmp added to support post-DFmax amplitude-change check.
void ParshlTracker::get_closest_frq(std::int64_t CurOsc,
                                   double DFmaxVal,
                                   std::int64_t Nlins,
                                   const std::vector<double>& LinFrq,
                                   const std::vector<double>& LinAmp) {
  // D6-observe quiet floor (precomputed once per call).
  // quiet_floor = 10^(-DAmax_observe/20): oscillators louder than this at a no-line stop
  // are counted as damax_exceeded_observe. Purely observational; no stop is added.
  const double damax_observe_floor =
      (p_.DAmax_observe > 0.0) ? std::pow(10.0, -p_.DAmax_observe / 20.0) : 0.0;

  struct RecGuard {
    RecGuard() {
      ++g_rec_depth;
      g_recursion_depth_max_frame = std::max<std::int64_t>(g_recursion_depth_max_frame, g_rec_depth);
      g_recursion_depth_max_acc = std::max<std::int64_t>(g_recursion_depth_max_acc, g_rec_depth);
    }
    ~RecGuard() { --g_rec_depth; }
  } guard;

  if (CurOsc < 1 || CurOsc > p_.MaxOscs) return;
  if (Nlins <= 0) {
    g_stop_reason = "no_best_line";
    ++g_no_best_line_frame;
    ++g_no_best_line_acc;
    stop_osc(CurOsc);
    return;
  }

  const double PF = st_.OscFrq[static_cast<std::size_t>(CurOsc)]; // previous target frequency

  if (trk_log_osc(CurOsc)) {
    std::ostringstream oss;
    oss << "match_begin frame1=" << g_frame1
        << " osc=" << CurOsc
        << " PF_Hz=" << PF
        << " DFmaxHz=" << DFmaxVal
        << " Nlins=" << Nlins
        << " recDepth=" << g_rec_depth;
    trk_log(oss.str());
  }

  std::int64_t BestLin = 0;
  double MinDist = std::numeric_limits<double>::infinity();

  for (std::int64_t CurLin = 1; CurLin <= Nlins; ++CurLin) {
    const double curFrq = get_lin(LinFrq, CurLin, Nlins);
    const double CurDist = absd(PF - curFrq);

    if (CurDist >= MinDist) continue;

    if (lin_free(CurLin)) {
      MinDist = CurDist;
      BestLin = CurLin;
    } else {
      // Claimed: resolve collision using PARSHL rule:
      // If I'm closer to the line than the current owner is, bump owner (recursive).
      const std::int64_t OtherOsc = st_.OscOfLin[static_cast<std::size_t>(CurLin)];
      if (OtherOsc < 1 || OtherOsc > p_.MaxOscs) continue;

      const double OtherDist = absd(curFrq - st_.OscFrq[static_cast<std::size_t>(OtherOsc)]);
      if (CurDist < OtherDist) {
        ++g_collision_count_frame;
        ++g_collision_count_acc;
        if (trk_verbose_enabled() && (trk_log_osc(CurOsc) || trk_log_osc(OtherOsc))) {
          std::ostringstream oss;
          oss << "collision frame1=" << g_frame1
              << " lin=" << CurLin
              << " curOsc=" << CurOsc
              << " otherOsc=" << OtherOsc
              << " curFrqHz=" << curFrq
              << " curDistHz=" << CurDist
              << " otherDistHz=" << OtherDist
              << " winner=curOsc"
              << " recDepth=" << g_rec_depth;
          trk_log(oss.str());
        }

        MinDist = CurDist;
        BestLin = CurLin;

        // Tentatively claim while the other guy fishes.
        const std::int64_t saved = st_.OscOfLin[static_cast<std::size_t>(CurLin)];
        const std::int64_t old_claim = st_.OscOfLin[static_cast<std::size_t>(CurLin)];
        st_.OscOfLin[static_cast<std::size_t>(CurLin)] = CurOsc;
        mark_oscoflin_write(CurLin, old_claim, CurOsc, "get_closest_frq", "collision_tentative_claim");

        // SAIL-faithful behavior: the displaced oscillator inherits the DFmax window of the
        // displacing oscillator during recursive reassignment.
        // See docs/TRACKER_DFMAX_RECURSION.md
        //
        // SAIL L1289: GetClosestFrq(OtherOsc, DFmax)
        //   — passes the CALLER's DFmax parameter, not dfmax(OtherOsc).
        //   The displaced oscillator must find a new line within the *displacer's*
        //   frequency-tolerance window, creating an asymmetry that depends on which
        //   oscillator initiated the chain:
        //     • Low-frequency displacer (narrow DFmax) → bumped osc has tight constraint → more stops.
        //     • High-frequency displacer (wide DFmax)  → bumped osc has generous constraint.
        //   This acts as a cascade brake: recursion terminates early when the
        //   displacer has a narrow window.  This is SAIL-faithful and intentionally preserved in V2.
        //   A possible V3 refinement would pass dfmax(OtherOsc) for symmetric behavior.
        get_closest_frq(OtherOsc, DFmaxVal, Nlins, LinFrq, LinAmp);

        // Unclaim (SAIL/PARSHL contract: set line owner to 0; final claim is applied in TakeIt)
        const std::int64_t old_restore = st_.OscOfLin[static_cast<std::size_t>(CurLin)];
        st_.OscOfLin[static_cast<std::size_t>(CurLin)] = 0;
        mark_oscoflin_write(CurLin, old_restore, 0, "get_closest_frq", "collision_restore_saved");
        ++g_collision_restore_saved_frame;
        ++g_collision_restore_saved_acc;

        const std::int64_t actual_restore = st_.OscOfLin[static_cast<std::size_t>(CurLin)];

        if (trk_checks_enabled()) {
          RestoreEvent ev;
          ev.lin = CurLin;
          ev.savedOsc = saved;
          ev.curOsc = CurOsc;
          ev.otherOsc = OtherOsc;
          g_restore_events_frame.push_back(ev);
        }

        if (trk_shadow_enabled() && actual_restore != 0) {
          ++g_shadow_unclaim_delta_frame;
          ++g_shadow_unclaim_delta_acc;
          if (trk_verbose_enabled() && (trk_log_osc(CurOsc) || trk_log_osc(OtherOsc))) {
            std::ostringstream sh;
            sh << "shadow_unclaim frame1=" << g_frame1
               << " lin=" << CurLin
               << " curOsc=" << CurOsc
               << " otherOsc=" << OtherOsc
              << " actual_restore=" << actual_restore
               << " sail_unclaim=0"
               << " delta=1";
            trk_log(sh.str());
          }
        }

        if (trk_verbose_enabled() && (trk_log_osc(CurOsc) || trk_log_osc(OtherOsc))) {
          std::ostringstream oss;
          oss << "collision_unclaim frame1=" << g_frame1
              << " lin=" << CurLin
              << " curOsc=" << CurOsc
              << " otherOsc=" << OtherOsc
              << " mode=unclaim_zero"
              << " restoredOsc=" << actual_restore;
          trk_log(oss.str());
        }
      }
    }
  }

  if (BestLin <= 0) {
    g_stop_reason = "no_best_line";
    ++g_no_best_line_frame;
    ++g_no_best_line_acc;
    // D6-observe (§25/JOS-faithful): no-line stop — oscillator lost all candidate lines.
    // JOS SAIL L1340: "the previous amplitude tells whether this is a reasonable assumption
    //   [that the line disappeared], but the information is not being put to use here."
    // Count as observed loud loss when oscillator amplitude exceeds the quiet floor.
    // Stop is unconditional; damax_exceeded_observe is purely observational.
    if (p_.DAmax_observe > 0.0 &&
        st_.PrvOscAmp[static_cast<std::size_t>(CurOsc)] > damax_observe_floor) {
      ++ts_.damax_exceeded_observe;
    }
    if (trk_log_osc(CurOsc)) {
      std::ostringstream oss;
      oss << "stop_check frame1=" << g_frame1
          << " osc=" << CurOsc
          << " bestLin=" << BestLin
          << " PF_Hz=" << PF
          << " BestFrqHz=0"
          << " MinDistHz=" << MinDist
          << " DFmaxHz=" << DFmaxVal
          << " decision=stop"
          << " reason=no_best_line";
      trk_log(oss.str());
    }
    stop_osc(CurOsc);
    return;
  }

  const double BestFrq = get_lin(LinFrq, BestLin, Nlins);
  const double ratio = (DFmaxVal > 0.0) ? (MinDist / DFmaxVal) : std::numeric_limits<double>::infinity();

  // SAIL L1305: IF MinDist > DFmax THEN StopOsc  (strict greater-than).
  // Equality (MinDist == DFmaxVal) falls through to TakeIt — the line IS accepted.
  // Effective acceptance criterion: MinDist <= DFmaxVal.
  if (MinDist > DFmaxVal) {
    // frequency tolerance exceeded -> squelch/stop
    ++g_dfmax_exceeded_frame;
    ++g_dfmax_exceeded_acc;
    ++ts_.dfmax_exceeded;
    g_stop_reason = "dfmax_exceeded";
    // D6-observe (§25/JOS-faithful): frequency-too-far is also a "no-line effective" stop.
    // Flag as observed loud loss when the oscillator was acoustically present.
    if (p_.DAmax_observe > 0.0 &&
        st_.PrvOscAmp[static_cast<std::size_t>(CurOsc)] > damax_observe_floor) {
      ++ts_.damax_exceeded_observe;
    }
    if (trk_log_osc(CurOsc)) {
      std::ostringstream oss;
      oss << "stop_check frame1=" << g_frame1
          << " osc=" << CurOsc
          << " bestLin=" << BestLin
          << " PF_Hz=" << PF
          << " BestFrqHz=" << BestFrq
          << " MinDistHz=" << MinDist
          << " DFmaxHz=" << DFmaxVal
          << " decision=stop"
          << " reason=dfmax_exceeded";
      trk_log(oss.str());
    }
    record_ratio(RatioRec{ratio, CurOsc, BestLin, PF, BestFrq, MinDist, DFmaxVal, "stop"});
    stop_osc(CurOsc);
    return;
  }

  // D6-observe (§25) checks are at the two stop paths above.
  //
  // D6-gate (§24/reconstructive/experimental): post-DFmax / pre-TakeIt amplitude gate.
  //
  // Archaeological citations:
  //   SAIL L158:  global declaration alongside DFmax.
  //   SAIL L962:  commented, never executed: "# DAmax <- 10; # Disabled;"
  //   JOS L1340 note: describes the no-line case, NOT this TakeIt case (see D6-observe).
  //
  // Reconstructive interpretation (extrapolated, NOT confirmed SAIL-executable):
  //   DAmax_gate = max tolerated amplitude change in dB on a frequency-matched line.
  //   Condition: |20*log10(LinAmp[BestLin] / PrvOscAmp[CurOsc])| > DAmax_gate => StopOsc.
  //   This ADDS stops beyond what frequency-matching alone would produce.
  //
  // Guard conditions:
  //   1. p_.DAmax_gate <= 0.0      -> feature disabled (default); skip entirely.
  //   2. PrvOscAmp[CurOsc] == 0   -> new/recycled/silent oscillator; bypass.
  //   3. LinAmp[BestLin] == 0     -> degenerate partial; bypass defensively.
  //
  // SigScl cancels in cur_amp/prv_amp; no access to SigScl needed.
  if (p_.DAmax_gate > 0.0) {
    const double prv_amp = st_.PrvOscAmp[static_cast<std::size_t>(CurOsc)];
    if (prv_amp > 0.0) {
      const double cur_amp = get_lin(LinAmp, BestLin, Nlins);
      if (cur_amp > 0.0) {
        const double amp_change_dB = std::abs(20.0 * std::log10(cur_amp / prv_amp));
        if (amp_change_dB > p_.DAmax_gate) {
          ++ts_.damax_gate_stops;
          g_stop_reason = "damax_gate";
          if (trk_log_osc(CurOsc)) {
            std::ostringstream oss;
            oss << "stop_check frame1=" << g_frame1
                << " osc=" << CurOsc
                << " bestLin=" << BestLin
                << " PF_Hz=" << PF
                << " BestFrqHz=" << BestFrq
                << " PrvAmpLin=" << prv_amp
                << " CurAmpLin=" << cur_amp
                << " AmpChangedB=" << amp_change_dB
                << " DAmax_gate=" << p_.DAmax_gate
                << " decision=stop"
                << " reason=damax_gate";
            trk_log(oss.str());
          }
          record_ratio(RatioRec{ratio, CurOsc, BestLin, PF, BestFrq, MinDist, DFmaxVal, "stop_damax_gate"});
          stop_osc(CurOsc);
          return;
        }
      }
    }
  }

  // TakeIt: assign
  const std::int64_t old_lin = st_.LinOfOsc[static_cast<std::size_t>(CurOsc)];
  const std::int64_t old_osc = st_.OscOfLin[static_cast<std::size_t>(BestLin)];
  st_.LinOfOsc[static_cast<std::size_t>(CurOsc)] = BestLin;
  st_.OscOfLin[static_cast<std::size_t>(BestLin)] = CurOsc;
  mark_linofosc_write(CurOsc, old_lin, BestLin, "get_closest_frq", "take_line_assign");
  mark_oscoflin_write(BestLin, old_osc, CurOsc, "get_closest_frq", "take_line_assign");
  if (CurOsc < static_cast<std::int64_t>(g_osc_action.size())) {
    g_osc_action[static_cast<std::size_t>(CurOsc)] = 1;
  }
  record_ratio(RatioRec{ratio, CurOsc, BestLin, PF, BestFrq, MinDist, DFmaxVal, "take"});

  if (trk_log_osc(CurOsc)) {
    std::ostringstream oss;
    oss << "take_line frame1=" << g_frame1
        << " osc=" << CurOsc
        << " PF_Hz=" << PF
        << " bestLin=" << BestLin
        << " BestFrqHz=" << BestFrq
        << " MinDistHz=" << MinDist
        << " DFmaxHz=" << DFmaxVal
        << " decision=take";
    trk_log(oss.str());
  }
}

// V2: Scan slots 1..Noscs for a dead oscillator that can be recycled.
// A slot is recyclable when it is not currently ON (LinOfOsc<=0 after Step 3 matching)
// AND both its current and previous amplitude outputs are zero (synthesis has fully decayed).
// Returns 1-based slot index, or 0 if none found.
std::int64_t ParshlTracker::find_free_osc() const noexcept {
  for (std::int64_t osc = 1; osc <= st_.Noscs; ++osc) {
    const std::size_t oi = static_cast<std::size_t>(osc);
    if (st_.LinOfOsc[oi] <= 0
        && st_.OscAmp[oi]    == 0.0
        && st_.PrvOscAmp[oi] == 0.0) {
      return osc;
    }
  }
  return 0;
}

// SAIL source: PROCEDURE UpdateMap
// Parshl-source.txt (tracker section, called from main analysis loop after FindPartials)
//
// Maps spectral partials (LinFrq/LinAmp) to synthesis oscillators each frame.
// Algorithm:
//   Frame1<=1 (first frame): direct assignment, no matching needed.
//   Normal frames — 5 steps:
//     Step 1: ARRTRAN(PrvLinOfOsc, LinOfOsc)  — copy current to previous
//     Step 2: ARRCLR(LinOfOsc), ARRCLR(OscOfLin)  — clear for this frame
//     Step 3: FOR CurOsc: GetClosestFrq(CurOsc, DFmax(CurOsc), Nlins, LinFrq)  (matching)
//     Step 4: FOR unclaimed Lin: assign to new oscillator  (add)
//     Step 5: OutAF  — update PrvOscFrq/PrvOscAmp/OscFrq/OscAmp for synthesis
//               SAIL L1548: IF PrvOscAmp[CurOsc]=0 THEN PrvOscFrq <- OscFrq  (FIX C-1)
void ParshlTracker::update_map(std::int64_t Frame1,
                              std::int64_t Nlins,
                              const std::vector<double>& LinAmp,
                              const std::vector<double>& LinFrq) {
  g_frame1 = Frame1;
  g_stop_reason = "unset";
  g_dfmax_exceeded_frame = 0;
  g_no_best_line_frame = 0;
  g_collision_count_frame = 0;
  g_recursion_depth_max_frame = 0;
  g_invariants_fail_frame = 0;
  g_stop_count_frame = 0;
  g_add_count_frame = 0;
  g_i2_fail_count_frame = 0;
  g_i11_fail_count_frame = 0;
  g_collision_restore_saved_frame = 0;
  g_shadow_unclaim_delta_frame = 0;
  g_shadow_instantrise_delta_frame = 0;
  g_dfmax_ratio_hist_frame = {{0, 0, 0, 0}};
  g_dfmax_top_frame.clear();
  g_restore_events_frame.clear();

  g_osc_action.assign(static_cast<std::size_t>(p_.MaxOscs + 1), 0);
  ensure_write_meta_sizes(p_.MaxOscs, p_.MaxLins);

  auto log_snapshot = [&](const char* phase) {
    if (!trk_snapshot_enabled() || !trk_frame_in_window()) return;

    std::int64_t on = 0, free = 0, squelched = 0, prv_nonzero = 0, claimed = 0;
    for (std::int64_t osc = 1; osc <= p_.MaxOscs; ++osc) {
      const std::int64_t lin = st_.LinOfOsc[static_cast<std::size_t>(osc)];
      if (lin > 0) ++on;
      else if (lin < 0) ++squelched;
      else ++free;
      if (st_.PrvLinOfOsc[static_cast<std::size_t>(osc)] != 0) ++prv_nonzero;
    }
    for (std::int64_t lin = 1; lin <= p_.MaxLins; ++lin) {
      if (st_.OscOfLin[static_cast<std::size_t>(lin)] > 0) ++claimed;
    }

    std::ostringstream oss;
    oss << "snapshot phase=" << phase
        << " frame1=" << Frame1
        << " Nlins=" << Nlins
        << " Noscs=" << st_.Noscs
        << " on=" << on
        << " free=" << free
        << " squelched=" << squelched
        << " claimed=" << claimed
        << " prv_nonzero=" << prv_nonzero
        << " hash_LinOfOsc=" << hash_i64_vec(st_.LinOfOsc)
        << " hash_PrvLinOfOsc=" << hash_i64_vec(st_.PrvLinOfOsc)
        << " hash_OscOfLin=" << hash_i64_vec(st_.OscOfLin)
        << " hash_OscFrq=" << hash_f64_vec(st_.OscFrq)
        << " hash_OscAmp=" << hash_f64_vec(st_.OscAmp)
        << " hash_PrvOscFrq=" << hash_f64_vec(st_.PrvOscFrq)
        << " hash_PrvOscAmp=" << hash_f64_vec(st_.PrvOscAmp);
    trk_log(oss.str());
  };

  // Defensive clamps: Nlins cannot exceed our storage.
  if (Nlins < 0) Nlins = 0;
  if (Nlins > p_.MaxLins) Nlins = p_.MaxLins;

  const std::int64_t Noscs_before = st_.Noscs;
  const std::int64_t want_pre = std::min<std::int64_t>(Nlins, p_.MaxOscs);
  const std::int64_t Nadd_pre = want_pre - Noscs_before;

  if (Frame1 <= 1 && trk_frame_in_window() && trk_verbose_enabled()) {
    std::int64_t nonzero_lin = 0;
    for (std::size_t i = 1; i < st_.LinOfOsc.size(); ++i) {
      if (st_.LinOfOsc[i] != 0) ++nonzero_lin;
    }
    std::int64_t nonzero_osc = 0;
    for (std::size_t i = 1; i < st_.OscOfLin.size(); ++i) {
      if (st_.OscOfLin[i] != 0) ++nonzero_osc;
    }
    std::ostringstream oss;
    oss << "init_probe frame1=" << Frame1
        << " Noscs=" << st_.Noscs
        << " LinOfOsc_size=" << st_.LinOfOsc.size()
        << " OscOfLin_size=" << st_.OscOfLin.size()
        << " LinOfOsc_nonzero_before_copy=" << nonzero_lin
        << " OscOfLin_nonzero_before_copy=" << nonzero_osc;
    trk_log(oss.str());
  }

  if (trk_frame_in_window()) {
    std::ostringstream oss;
    oss << "frame_begin frame1=" << Frame1
        << " frame0=" << (Frame1 - 1)
        << " Nlins=" << Nlins
        << " Noscs_before=" << Noscs_before
        << " want=" << want_pre
        << " Nadd=" << Nadd_pre
        << " InstantRise=" << (p_.InstantRise ? 1 : 0)
        << " MaxOscs=" << p_.MaxOscs
        << " Fc1=" << p_.Fc1
        << " Fc2=" << p_.Fc2
        << " DFmax1=" << p_.DFmax1
        << " DFmax2=" << p_.DFmax2;
    trk_log(oss.str());
  }

  log_snapshot("begin");

  // SAIL: Frame1<=1 — first-frame initialisation (direct assignment, no matching).
  // All oscillators are allocated sequentially to the Nlins found partials.
  // SAIL L1481: PrvOscAmp[CurOsc] <- (IF InstantRise THEN LinAmp[CurLin] ELSE 0)
  //   InstantRise=FALSE (default): PrvOscAmp starts at 0, synthesiser ramps from silence.
  //   InstantRise=TRUE: PrvOscAmp = LinAmp, synthesiser starts at full amplitude.
  if (Frame1 <= 1) {
    st_.Noscs = std::min<std::int64_t>(Nlins, p_.MaxOscs);

    if (trk_trace_writes_enabled()) {
      for (std::int64_t osc = 1; osc <= p_.MaxOscs; ++osc) {
        const std::size_t oi = static_cast<std::size_t>(osc);
        const std::int64_t old_lin = st_.LinOfOsc[oi];
        st_.LinOfOsc[oi] = 0;
        mark_linofosc_write(osc, old_lin, 0, "update_map", "init_clear");
        st_.PrvLinOfOsc[oi] = 0;
      }
      for (std::int64_t lin = 1; lin <= p_.MaxLins; ++lin) {
        const std::size_t li = static_cast<std::size_t>(lin);
        const std::int64_t old_osc = st_.OscOfLin[li];
        st_.OscOfLin[li] = 0;
        mark_oscoflin_write(lin, old_osc, 0, "update_map", "init_clear");
      }
    } else {
      std::fill(st_.LinOfOsc.begin(), st_.LinOfOsc.end(), 0);
      std::fill(st_.PrvLinOfOsc.begin(), st_.PrvLinOfOsc.end(), 0);
      std::fill(st_.OscOfLin.begin(), st_.OscOfLin.end(), 0);
    }

    for (std::int64_t osc = 1; osc <= st_.Noscs; ++osc) {
      const std::size_t oi = static_cast<std::size_t>(osc);
      const std::int64_t old_lin = st_.LinOfOsc[oi];
      const std::int64_t old_osc = st_.OscOfLin[oi];
      st_.LinOfOsc[oi] = osc;
      st_.OscOfLin[oi] = osc;
      st_.OscFrq[oi] = get_lin(LinFrq, osc, Nlins);
      st_.OscAmp[oi] = get_lin(LinAmp, osc, Nlins);
      st_.PrvOscFrq[oi] = get_lin(LinFrq, osc, Nlins);
      // FIX A: SAIL line 1481: PrvOscAmp[CurOsc] <- (IF InstantRise THEN LinAmp[CurLin] ELSE 0)
      // Default InstantRise=FALSE → PrvOscAmp starts at 0, synth ramps from silence.
      st_.PrvOscAmp[oi] = p_.InstantRise ? get_lin(LinAmp, osc, Nlins) : 0.0;
      if (trk_trace_writes_enabled()) {
        mark_linofosc_write(osc, old_lin, osc, "update_map", "init_assign");
        mark_oscoflin_write(osc, old_osc, osc, "update_map", "init_assign");
      }
    }

    ts_.add_count += st_.Noscs;  // first-frame init: all Noscs oscillators are new starts

    if (trk_frame_in_window()) {
      std::ostringstream oss;
      oss << "init_frame_done frame1=" << Frame1
          << " Nlins=" << Nlins
          << " Noscs=" << st_.Noscs;
      trk_log(oss.str());
    }

    log_snapshot("end");
    return;
  }

  // Step 1: SAIL ARRTRAN(PrvLinOfOsc, LinOfOsc) — copy current mapping to previous
  st_.PrvLinOfOsc = st_.LinOfOsc;

  // Step 2: SAIL ARRCLR(LinOfOsc), ARRCLR(OscOfLin) — clear maps for new frame
  if (trk_trace_writes_enabled()) {
    for (std::int64_t osc = 1; osc <= p_.MaxOscs; ++osc) {
      const std::size_t oi = static_cast<std::size_t>(osc);
      const std::int64_t oldv = st_.LinOfOsc[oi];
      st_.LinOfOsc[oi] = 0;
      mark_linofosc_write(osc, oldv, 0, "update_map", "clear_maps");
    }
    for (std::int64_t lin = 1; lin <= p_.MaxLins; ++lin) {
      const std::size_t li = static_cast<std::size_t>(lin);
      const std::int64_t oldv = st_.OscOfLin[li];
      st_.OscOfLin[li] = 0;
      mark_oscoflin_write(lin, oldv, 0, "update_map", "clear_maps");
    }
  } else {
    std::fill(st_.LinOfOsc.begin(), st_.LinOfOsc.end(), 0);
    std::fill(st_.OscOfLin.begin(), st_.OscOfLin.end(), 0);
  }

  if (trk_frame_in_window()) {
    const std::int64_t osc_lim = trk_osc_max();
    const std::int64_t ko = (osc_lim <= 0) ? p_.MaxOscs : std::min<std::int64_t>(osc_lim, p_.MaxOscs);
    trk_log(std::string("maps_after_copy_clear prvLin=") + dump_i64_1based(st_.PrvLinOfOsc, ko));
    trk_log(std::string("maps_after_copy_clear lin=") + dump_i64_1based(st_.LinOfOsc, ko));
  }

  // I4: post-clear integrity (must be zero before matching/add)
  {
    std::int64_t bad_lin = 0;
    for (std::size_t i = 1; i < st_.LinOfOsc.size(); ++i) {
      if (st_.LinOfOsc[i] != 0) { ++bad_lin; }
    }
    std::int64_t bad_osc = 0;
    for (std::size_t i = 1; i < st_.OscOfLin.size(); ++i) {
      if (st_.OscOfLin[i] != 0) { ++bad_osc; }
    }
    if (bad_lin != 0 || bad_osc != 0) {
      ++g_invariants_fail_frame;
      ++g_invariants_fail_acc;
      if (trk_frame_in_window()) {
        std::ostringstream oss;
        oss << "invariant I4 FAIL frame1=" << Frame1
            << " bad_LinOfOsc_nonzero=" << bad_lin
            << " bad_OscOfLin_nonzero=" << bad_osc;
        trk_log(oss.str());
      }
    } else if (trk_frame_in_window() && trk_verbose_enabled()) {
      trk_log("invariant I4 PASS (post-clear maps are zero)");
    }
  }

  // I7: monotonicity of LinFrq for active lines
  {
    if (Nlins <= 1) {
      if (trk_frame_in_window() && trk_verbose_enabled()) {
        trk_log("invariant I7 SKIP (Nlins<=1)");
      }
    } else {
      bool ordered = true;
      std::int64_t first_bad = -1;
      double prev = get_lin(LinFrq, 1, Nlins);
      for (std::int64_t lin = 2; lin <= Nlins; ++lin) {
        const double cur = get_lin(LinFrq, lin, Nlins);
        if (cur < prev) {
          ordered = false;
          first_bad = lin;
          break;
        }
        prev = cur;
      }
      if (!ordered) {
        ++g_invariants_fail_frame;
        ++g_invariants_fail_acc;
        if (trk_frame_in_window()) {
          std::ostringstream oss;
          oss << "invariant I7 FAIL frame1=" << Frame1
              << " first_out_of_order_lin=" << first_bad
              << " prevLinFrqHz=" << get_lin(LinFrq, first_bad - 1, Nlins)
              << " curLinFrqHz=" << get_lin(LinFrq, first_bad, Nlins);
          trk_log(oss.str());
        }
      } else if (trk_frame_in_window() && trk_verbose_enabled()) {
        trk_log("invariant I7 PASS (LinFrq nondecreasing)");
      }
    }
  }

  if (trk_frame_in_window() && trk_verbose_enabled()) {
    const std::int64_t lim = trk_lin_max();
    const std::int64_t k = (lim <= 0) ? Nlins : std::min<std::int64_t>(lim, Nlins);
    for (std::int64_t lin = 1; lin <= k; ++lin) {
      std::ostringstream oss;
      oss << "lin_preview frame1=" << Frame1
          << " lin=" << lin
          << " frqHz=" << get_lin(LinFrq, lin, Nlins)
          << " amp=" << get_lin(LinAmp, lin, Nlins);
      trk_log(oss.str());
    }
  }

  // If no lines, stop all allocated oscillators and update targets accordingly.
  if (Nlins == 0) {
    for (std::int64_t osc = 1; osc <= st_.Noscs; ++osc) stop_osc(osc);
  } else {
    // Step 3: SAIL matching — FOR CurOsc: GetClosestFrq(CurOsc, DFmax(CurOsc), Nlins, LinFrq)
    // Assigns each existing oscillator to the nearest unclaimed partial within DFmax.
    // D6: LinAmp forwarded for the gate variant; also used by observe precomputation guard.
    for (std::int64_t CurOsc = 1; CurOsc <= st_.Noscs; ++CurOsc) {
      get_closest_frq(CurOsc, dfmax(CurOsc), Nlins, LinFrq, LinAmp);
    }

    // Step 4: SAIL Add — allocate new oscillators to unclaimed spectral lines
    // FOR unclaimed CurLin: OscOfLin[CurLin]=CurOsc, LinOfOsc[CurOsc]=CurLin
    // Legacy mode: exact SAIL behaviour — Noscs is monotone non-decreasing.
    // V2 mode:     try find_free_osc() (dead silent slot) before extending Noscs;
    //              this prevents bank saturation on dense signals (D5).
    const std::int64_t want = std::min<std::int64_t>(Nlins, p_.MaxOscs);
    const std::int64_t Nadd = want - st_.Noscs;

    if (trk_frame_in_window()) {
      std::ostringstream oss;
      oss << "add_phase frame1=" << Frame1
          << " want=" << want
          << " Noscs_before_add=" << st_.Noscs
          << " Nadd=" << Nadd;
      trk_log(oss.str());
    }

    if (p_.legacy_mode) {
      // ── LEGACY: exact SAIL monotone-Noscs, extend-only ───────────────────
      if (Nadd > 0) {
        std::int64_t CurLin = 0;
        const std::int64_t oldNoscs = st_.Noscs;
        for (std::int64_t CurOsc = st_.Noscs + 1; CurOsc <= st_.Noscs + Nadd; ++CurOsc) {
          do {
            ++CurLin;
            if (CurLin > Nlins) break;
          } while (st_.OscOfLin[static_cast<std::size_t>(CurLin)] != 0);

          if (CurLin > Nlins) break;

          const std::int64_t old_osc = st_.OscOfLin[static_cast<std::size_t>(CurLin)];
          const std::int64_t old_lin = st_.LinOfOsc[static_cast<std::size_t>(CurOsc)];
          st_.OscOfLin[static_cast<std::size_t>(CurLin)] = CurOsc;
          st_.LinOfOsc[static_cast<std::size_t>(CurOsc)] = CurLin;
          mark_oscoflin_write(CurLin, old_osc, CurOsc, "update_map", "add_assign");
          mark_linofosc_write(CurOsc, old_lin, CurLin, "update_map", "add_assign");
          ++g_add_count_frame;
          ++g_add_count_acc;
          ++ts_.add_count;
          if (CurOsc < static_cast<std::int64_t>(g_osc_action.size())) {
            g_osc_action[static_cast<std::size_t>(CurOsc)] = 3;
          }

          if (trk_log_osc(CurOsc)) {
            std::ostringstream oss;
            oss << "add_assign frame1=" << Frame1
                << " osc=" << CurOsc
                << " lin=" << CurLin
                << " linFrqHz=" << get_lin(LinFrq, CurLin, Nlins)
                << " linAmp=" << get_lin(LinAmp, CurLin, Nlins);
            trk_log(oss.str());
          }
        }
        st_.Noscs = std::min<std::int64_t>(p_.MaxOscs, st_.Noscs + Nadd);

        if (trk_frame_in_window()) {
          std::ostringstream oss;
          oss << "add_done frame1=" << Frame1
              << " Noscs_before=" << oldNoscs
              << " Noscs_after=" << st_.Noscs;
          trk_log(oss.str());
        }
      }
    } else {
      // ── V2: recycle dead slots before extending Noscs ────────────────────
      // Iterate all unclaimed lines; for each, recycle a dead slot or extend.
      const std::int64_t oldNoscs = st_.Noscs;
      std::int64_t CurLin = 0;
      std::int64_t recycle_count = 0;

      while (true) {
        // Advance to next unclaimed line.
        bool found_lin = false;
        do {
          ++CurLin;
          if (CurLin > Nlins) break;
          if (st_.OscOfLin[static_cast<std::size_t>(CurLin)] == 0) { found_lin = true; break; }
        } while (true);
        if (!found_lin) break;

        // Get a slot: prefer recycling a dead one, then extend if room.
        std::int64_t CurOsc;
        bool recycled = false;

        const std::int64_t free_slot = find_free_osc();
        if (free_slot > 0) {
          CurOsc   = free_slot;
          recycled = true;
          ++recycle_count;
        } else if (st_.Noscs < p_.MaxOscs) {
          ++st_.Noscs;
          CurOsc = st_.Noscs;
        } else {
          // Fully saturated — no dead slots and no room to extend.
          break;
        }

        // Assign CurLin ↔ CurOsc.
        const std::int64_t old_osc = st_.OscOfLin[static_cast<std::size_t>(CurLin)];
        const std::int64_t old_lin = st_.LinOfOsc[static_cast<std::size_t>(CurOsc)];
        st_.OscOfLin[static_cast<std::size_t>(CurLin)] = CurOsc;
        st_.LinOfOsc[static_cast<std::size_t>(CurOsc)] = CurLin;
        mark_oscoflin_write(CurLin, old_osc, CurOsc, "update_map",
                            recycled ? "recycle_assign" : "add_assign");
        mark_linofosc_write(CurOsc, old_lin, CurLin, "update_map",
                            recycled ? "recycle_assign" : "add_assign");
        ++g_add_count_frame;
        ++g_add_count_acc;
        if (recycled) ++ts_.recycle_count; else ++ts_.add_count;
        if (CurOsc < static_cast<std::int64_t>(g_osc_action.size())) {
          g_osc_action[static_cast<std::size_t>(CurOsc)] = 3;
        }

        if (trk_log_osc(CurOsc)) {
          std::ostringstream oss;
          oss << (recycled ? "recycle_assign" : "add_assign")
              << " frame1=" << Frame1
              << " osc=" << CurOsc
              << " lin=" << CurLin
              << " linFrqHz=" << get_lin(LinFrq, CurLin, Nlins)
              << " linAmp=" << get_lin(LinAmp, CurLin, Nlins)
              << " recycled=" << (recycled ? 1 : 0);
          trk_log(oss.str());
        }
      }

      if (trk_frame_in_window()) {
        std::ostringstream oss;
        oss << "add_done frame1=" << Frame1
            << " Noscs_before=" << oldNoscs
            << " Noscs_after=" << st_.Noscs
            << " recycles=" << recycle_count;
        trk_log(oss.str());
      }
    }
  }

  // 5) Update targets for synthesis up to latest spectral frame (PARSHL OutAF order)
  for (std::int64_t CurOsc = 1; CurOsc <= p_.MaxOscs; ++CurOsc) {
    const std::size_t i = static_cast<std::size_t>(CurOsc);

    // Step 1: previous target frequency
    st_.PrvOscFrq[i] = st_.OscFrq[i];

    // Step 2: new target frequency if ON
    if (osc_on(CurOsc)) {
      const std::int64_t lin = st_.LinOfOsc[i];
      st_.OscFrq[i] = get_lin(LinFrq, lin, Nlins);
    }

    // Step 3: previous amplitude
    st_.PrvOscAmp[i] = st_.OscAmp[i];

    // Step 4: new amplitude
    if (osc_on(CurOsc)) {
      const std::int64_t lin = st_.LinOfOsc[i];
      st_.OscAmp[i] = get_lin(LinAmp, lin, Nlins);
    } else {
      st_.OscAmp[i] = 0.0;
    }

    // Step 5: SAIL OutAF (Parshl-source.txt L1548) — UNCONDITIONAL rule:
    //   IF PrvOscAmp[CurOsc]=0 THEN PrvOscFrq[CurOsc] <- OscFrq[CurOsc]
    // InstantRise does NOT appear in OutAF in SAIL; only in the Frame 1 initialisation
    // (L1481), a completely different context. The previous guard on
    // p_.InstantRise was factually incorrect (audit C-1).
    bool instant_rise_applied = false;
    if (st_.PrvOscAmp[i] == 0.0) {
      st_.PrvOscFrq[i] = st_.OscFrq[i];
      instant_rise_applied = true;
    }

    if (trk_log_osc(CurOsc)) {
      std::ostringstream oss;
      oss << "targets frame1=" << Frame1
          << " osc=" << CurOsc
          << " lin=" << st_.LinOfOsc[i]
          << " prvAmp=" << st_.PrvOscAmp[i]
          << " amp=" << st_.OscAmp[i]
          << " prvFrqHz=" << st_.PrvOscFrq[i]
          << " frqHz=" << st_.OscFrq[i]
          << " instantRiseApplied=" << (instant_rise_applied ? 1 : 0);
      trk_log(oss.str());
    }
  }

  // I1/I2/I3/I5/I6 checks (non-fatal)
  {
    std::int64_t fail_i1 = 0;
    std::int64_t fail_i2 = 0;
    std::int64_t fail_i3 = 0;
    std::int64_t fail_i5 = 0;
    std::int64_t fail_i6 = 0;
    std::int64_t fail_i8 = 0;
    std::int64_t fail_i9 = 0;
    std::int64_t fail_i9_prev0 = 0;
    std::int64_t fail_i10 = 0;
    std::int64_t fail_i11 = 0;
    std::int64_t i2_cases_logged = 0;

    const auto lin_class = [](std::int64_t lin) -> const char* {
      if (lin > 0) return "on";
      if (lin < 0) return "squelched";
      return "free";
    };

    std::vector<std::int64_t> seen_lin(static_cast<std::size_t>(std::max<std::int64_t>(1, Nlins + 1)), 0);
    std::vector<std::int64_t> seen_osc(static_cast<std::size_t>(std::max<std::int64_t>(1, st_.Noscs + 1)), 0);

    for (std::int64_t osc = 1; osc <= st_.Noscs; ++osc) {
      const std::size_t oi = static_cast<std::size_t>(osc);
      const std::int64_t lin = st_.LinOfOsc[oi];

      if (lin > p_.MaxLins || lin < -p_.MaxLins) ++fail_i6;
      if (std::llabs(lin) > Nlins) ++fail_i10;

      if (lin > 0) {
        if (lin > Nlins) ++fail_i10;
        if (lin >= static_cast<std::int64_t>(st_.OscOfLin.size()) ||
            st_.OscOfLin[static_cast<std::size_t>(lin)] != osc) {
          ++fail_i1;
        }
        if (lin >= 1 && lin < static_cast<std::int64_t>(seen_lin.size())) {
          if (seen_lin[static_cast<std::size_t>(lin)] != 0) ++fail_i8;
          seen_lin[static_cast<std::size_t>(lin)] = osc;
        }
      }

      if (lin < 0) {
        const std::int64_t prv = st_.PrvLinOfOsc[oi];
        if (std::llabs(lin) != std::llabs(prv)) ++fail_i3;
        if (prv == 0) {
          ++fail_i9_prev0;
          if (trk_frame_in_window() && fail_i9_prev0 <= trk_case_max()) {
            std::ostringstream i9;
            i9 << "I9_case frame1=" << Frame1
               << " osc=" << osc
               << " lin=" << lin
               << " prvLin=" << prv
               << " note=squelched_with_prev0";
            trk_log(i9.str());
          }
        }
        const std::int64_t abs_lin = std::llabs(lin);
        if (abs_lin >= 1 && abs_lin < static_cast<std::int64_t>(st_.OscOfLin.size()) &&
            st_.OscOfLin[static_cast<std::size_t>(abs_lin)] == osc) {
          ++fail_i9;
        }
      }

      if (!osc_on(osc) && st_.OscAmp[oi] != 0.0) ++fail_i5;
    }

    for (std::int64_t lin = 1; lin <= Nlins; ++lin) {
      const std::size_t li = static_cast<std::size_t>(lin);
      const std::int64_t osc = st_.OscOfLin[li];
      if (osc > p_.MaxOscs || osc < 0) ++fail_i6;
      if (osc > st_.Noscs) ++fail_i10;
      if (osc > 0) {
        if (osc >= 1 && osc < static_cast<std::int64_t>(seen_osc.size())) {
          if (seen_osc[static_cast<std::size_t>(osc)] != 0) ++fail_i8;
          seen_osc[static_cast<std::size_t>(osc)] = lin;
        }
        if (osc >= static_cast<std::int64_t>(st_.LinOfOsc.size()) ||
            st_.LinOfOsc[static_cast<std::size_t>(osc)] != lin) {
          ++fail_i2;
          if (trk_frame_in_window() && i2_cases_logged < trk_case_max()) {
            const std::int64_t lin_of_osc =
              (osc >= 1 && osc < static_cast<std::int64_t>(st_.LinOfOsc.size()))
                ? st_.LinOfOsc[static_cast<std::size_t>(osc)]
                : std::numeric_limits<std::int64_t>::min();
            const std::int64_t prv_lin_of_osc =
              (osc >= 1 && osc < static_cast<std::int64_t>(st_.PrvLinOfOsc.size()))
                ? st_.PrvLinOfOsc[static_cast<std::size_t>(osc)]
                : std::numeric_limits<std::int64_t>::min();

            const char* action = "none";
            if (osc >= 1 && osc < static_cast<std::int64_t>(g_osc_action.size())) {
              const int a = g_osc_action[static_cast<std::size_t>(osc)];
              action = (a == 1) ? "take" : (a == 2) ? "stop" : (a == 3) ? "add" : "none";
            }

            std::string w_lin = "{unavailable}";
            if (osc >= 1 && osc < static_cast<std::int64_t>(g_lastwrite_linofosc.size())) {
              w_lin = fmt_write_meta(g_lastwrite_linofosc[static_cast<std::size_t>(osc)]);
            }
            std::string w_osc = "{unavailable}";
            if (lin >= 1 && lin < static_cast<std::int64_t>(g_lastwrite_oscoflin.size())) {
              w_osc = fmt_write_meta(g_lastwrite_oscoflin[static_cast<std::size_t>(lin)]);
            }

            std::ostringstream c;
            c << "I2_case frame1=" << Frame1
              << " lin=" << lin
              << " osc=" << osc
              << " OscOfLin[lin]=" << osc
              << " LinOfOsc[osc]=" << lin_of_osc
              << " PrvLinOfOsc[osc]=" << prv_lin_of_osc
              << " endAction=" << action
              << " lastWrite_LinOfOsc=" << w_lin
              << " lastWrite_OscOfLin=" << w_osc;
            trk_log(c.str());

            const std::int64_t final_lin_of_osc =
              (osc >= 1 && osc < static_cast<std::int64_t>(st_.LinOfOsc.size()))
                ? st_.LinOfOsc[static_cast<std::size_t>(osc)]
                : std::numeric_limits<std::int64_t>::min();
            const std::int64_t final_osc_of_lin =
              (lin >= 1 && lin < static_cast<std::int64_t>(st_.OscOfLin.size()))
                ? st_.OscOfLin[static_cast<std::size_t>(lin)]
                : std::numeric_limits<std::int64_t>::min();

            std::ostringstream p;
            p << "I2_post frame1=" << Frame1
              << " phase=post_frame_final_after_update_outAF"
              << " lin=" << lin
              << " osc=" << osc
              << " LinOfOsc_final=" << final_lin_of_osc
              << " LinOfOsc_class=" << lin_class(final_lin_of_osc)
              << " OscOfLin_final=" << final_osc_of_lin;
            if (final_osc_of_lin > 0 && final_osc_of_lin != osc &&
                final_osc_of_lin < static_cast<std::int64_t>(st_.LinOfOsc.size())) {
              p << " OscOfLin_points_to_other_osc=" << final_osc_of_lin
                << " other_LinOfOsc=" << st_.LinOfOsc[static_cast<std::size_t>(final_osc_of_lin)];
            }
            trk_log(p.str());
            ++i2_cases_logged;
          }
        }
      }
    }

    for (const auto& ev : g_restore_events_frame) {
      if (ev.savedOsc <= 0 || ev.lin <= 0 || ev.lin > Nlins) continue;
      if (ev.savedOsc >= static_cast<std::int64_t>(g_osc_action.size())) continue;
      if (ev.lin >= static_cast<std::int64_t>(st_.OscOfLin.size())) continue;
      const bool saved_stopped = (g_osc_action[static_cast<std::size_t>(ev.savedOsc)] == 2);
      const bool saved_took_other = (g_osc_action[static_cast<std::size_t>(ev.savedOsc)] == 1) &&
        (ev.savedOsc < static_cast<std::int64_t>(st_.LinOfOsc.size())) &&
        (st_.LinOfOsc[static_cast<std::size_t>(ev.savedOsc)] != ev.lin);
      const bool line_still_claims_saved =
        (st_.OscOfLin[static_cast<std::size_t>(ev.lin)] == ev.savedOsc);
      if ((saved_stopped && line_still_claims_saved) || (saved_took_other && line_still_claims_saved)) {
        ++fail_i11;
        if (trk_frame_in_window() && fail_i11 <= trk_case_max()) {
          std::ostringstream i11;
          i11 << "I11_case frame1=" << Frame1
              << " lin=" << ev.lin
              << " savedOsc=" << ev.savedOsc
              << " curOsc=" << ev.curOsc
              << " otherOsc=" << ev.otherOsc
                << " savedEndAction=" << (saved_stopped ? "stop" : "take_other_lin")
                << " savedLin="
                << ((ev.savedOsc >= 1 && ev.savedOsc < static_cast<std::int64_t>(st_.LinOfOsc.size()))
                  ? st_.LinOfOsc[static_cast<std::size_t>(ev.savedOsc)]
                  : std::numeric_limits<std::int64_t>::min())
              << " OscOfLin[lin]=" << st_.OscOfLin[static_cast<std::size_t>(ev.lin)];
          trk_log(i11.str());

          const std::int64_t post_osc = st_.OscOfLin[static_cast<std::size_t>(ev.lin)];
          const std::int64_t post_lin_saved =
            (ev.savedOsc >= 1 && ev.savedOsc < static_cast<std::int64_t>(st_.LinOfOsc.size()))
              ? st_.LinOfOsc[static_cast<std::size_t>(ev.savedOsc)]
              : std::numeric_limits<std::int64_t>::min();
          std::ostringstream p11;
          p11 << "I11_post frame1=" << Frame1
              << " phase=post_frame_final_after_update_outAF"
              << " lin=" << ev.lin
              << " osc=" << ev.savedOsc
              << " LinOfOsc_final=" << post_lin_saved
              << " LinOfOsc_class=" << lin_class(post_lin_saved)
              << " OscOfLin_final=" << post_osc;
          if (post_osc > 0 && post_osc != ev.savedOsc &&
              post_osc < static_cast<std::int64_t>(st_.LinOfOsc.size())) {
            p11 << " OscOfLin_points_to_other_osc=" << post_osc
                << " other_LinOfOsc=" << st_.LinOfOsc[static_cast<std::size_t>(post_osc)];
          }
          trk_log(p11.str());
        }
      }
    }

    g_i2_fail_count_frame = fail_i2;
    g_i11_fail_count_frame = fail_i11;
    g_i2_fail_count_acc += fail_i2;
    g_i11_fail_count_acc += fail_i11;

    const std::int64_t fail_sum =
      fail_i1 + fail_i2 + fail_i3 + fail_i5 + fail_i6 + fail_i8 + fail_i9 + fail_i10 + fail_i11;
    g_invariants_fail_frame += fail_sum;
    g_invariants_fail_acc += fail_sum;

    if (trk_frame_in_window()) {
      std::ostringstream oss;
      oss << "invariants frame1=" << Frame1
          << " I1=" << (fail_i1 == 0 ? "PASS" : "FAIL(" + std::to_string(fail_i1) + ")")
          << " I2=" << (fail_i2 == 0 ? "PASS" : "FAIL(" + std::to_string(fail_i2) + ")")
          << " I3=" << (fail_i3 == 0 ? "PASS" : "FAIL(" + std::to_string(fail_i3) + ")")
          << " I5=" << (fail_i5 == 0 ? "PASS" : "FAIL(" + std::to_string(fail_i5) + ")")
          << " I6=" << (fail_i6 == 0 ? "PASS" : "FAIL(" + std::to_string(fail_i6) + ")")
          << " I8=" << (fail_i8 == 0 ? "PASS" : "FAIL(" + std::to_string(fail_i8) + ")")
          << " I9=" << ((fail_i9 == 0 && fail_i9_prev0 == 0)
                         ? "PASS"
                         : "FAIL(" + std::to_string(fail_i9 + fail_i9_prev0) + ")")
          << " I10=" << (fail_i10 == 0 ? "PASS" : "FAIL(" + std::to_string(fail_i10) + ")")
          << " I11=" << (fail_i11 == 0 ? "PASS" : "FAIL(" + std::to_string(fail_i11) + ")")
          << " frame_fail_total=" << g_invariants_fail_frame;
      trk_log(oss.str());
    }
  }

  log_snapshot("end");

  if (trk_frame_in_window()) {
    std::ostringstream oss;
    oss << "counters frame1=" << Frame1
        << " dfmax_exceeded_frame=" << g_dfmax_exceeded_frame
        << " no_best_line_frame=" << g_no_best_line_frame
      << " i2_fail_count_frame=" << g_i2_fail_count_frame
      << " i11_fail_count_frame=" << g_i11_fail_count_frame
      << " collision_restore_saved_frame=" << g_collision_restore_saved_frame
        << " stop_count_frame=" << g_stop_count_frame
        << " add_count_frame=" << g_add_count_frame
        << " collision_count_frame=" << g_collision_count_frame
        << " recursion_depth_max_frame=" << g_recursion_depth_max_frame
        << " shadow_unclaim_delta_frame=" << g_shadow_unclaim_delta_frame
        << " shadow_instantrise_delta_frame=" << g_shadow_instantrise_delta_frame
        << " invariants_fail_frame=" << g_invariants_fail_frame
        << " | accum dfmax_exceeded=" << g_dfmax_exceeded_acc
        << " no_best_line=" << g_no_best_line_acc
        << " i2_fail_count=" << g_i2_fail_count_acc
        << " i11_fail_count=" << g_i11_fail_count_acc
        << " collision_restore_saved=" << g_collision_restore_saved_acc
        << " stop_count=" << g_stop_count_acc
        << " add_count=" << g_add_count_acc
        << " collision_count=" << g_collision_count_acc
        << " recursion_depth_max=" << g_recursion_depth_max_acc
        << " shadow_unclaim_delta=" << g_shadow_unclaim_delta_acc
        << " shadow_instantrise_delta=" << g_shadow_instantrise_delta_acc
        << " invariants_fail_count=" << g_invariants_fail_acc;
    trk_log(oss.str());

    if (trk_dfmax_stats_enabled()) {
      std::ostringstream hs;
      hs << "dfmax_ratio_hist frame1=" << Frame1
         << " le1=" << g_dfmax_ratio_hist_frame[0]
         << " b1_2=" << g_dfmax_ratio_hist_frame[1]
         << " b2_5=" << g_dfmax_ratio_hist_frame[2]
         << " gt5=" << g_dfmax_ratio_hist_frame[3]
         << " | accum le1=" << g_dfmax_ratio_hist_acc[0]
         << " b1_2=" << g_dfmax_ratio_hist_acc[1]
         << " b2_5=" << g_dfmax_ratio_hist_acc[2]
         << " gt5=" << g_dfmax_ratio_hist_acc[3];
      trk_log(hs.str());

      for (std::size_t i = 0; i < g_dfmax_top_frame.size(); ++i) {
        const auto& r = g_dfmax_top_frame[i];
        std::ostringstream tk;
        tk << "dfmax_ratio_top frame1=" << Frame1
           << " rank=" << (i + 1)
           << " ratio=" << r.ratio
           << " osc=" << r.osc
           << " bestLin=" << r.bestLin
           << " PF_Hz=" << r.PF
           << " BestFrqHz=" << r.BestFrq
           << " MinDistHz=" << r.MinDist
           << " DFmaxHz=" << r.DFmax
           << " decision=" << r.decision;
        trk_log(tk.str());
      }
    }

    if (trk_shadow_enabled()) {
      std::ostringstream sh;
      sh << "shadow_summary frame1=" << Frame1
         << " unclaim_delta_frame=" << g_shadow_unclaim_delta_frame
         << " instantrise_delta_frame=" << g_shadow_instantrise_delta_frame
         << " dfmax_exceeded_frame=" << g_dfmax_exceeded_frame
         << " invariants_fail_frame=" << g_invariants_fail_frame;
      trk_log(sh.str());
    }
  }
}

} // namespace parshl
