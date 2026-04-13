// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// parshl_peak_claim.cpp
// [V4 T1] Peak-claiming tracker — minimal prototype.
//
// Origin:  E5 failure analysis (§47C); Architecture A1 (docs/PARSHL_V3_FINAL.md §8).
// Design:  docs/T1_PEAK_CLAIMING_TRACKER.md (§54).
// Impl:    §55.
//
// See parshl_peak_claim.hpp for the public API and contract.

#include "parshl_peak_claim.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

namespace parshl {

T1Stats t1_update(
    const TrackerState&          prv,
    const TrackerParams&         tp,
    std::int64_t                 /*Frame1*/,
    std::int64_t                 Nlins,
    const std::vector<double>&   LinAmp,
    const std::vector<double>&   LinFrq,
    int                          squelch_hold,
    std::vector<int>&            hold_counts,
    double                       t1_amp_hysteresis,  // [§64]
    TrackerState&                out)
{
  T1Stats stats{};

  const auto MaxOscs = tp.MaxOscs;
  const std::size_t sz  = static_cast<std::size_t>(MaxOscs + 1);  // 1-based size
  const std::size_t lsz = static_cast<std::size_t>(tp.MaxLins + 1);

  // ── 1. Initialise out: copy Prv* from prv, zero current fields ───────────
  out.Noscs = 0;
  out.df_ready = prv.df_ready;
  out.alpha    = prv.alpha;
  out.beta     = prv.beta;

  // Resize and zero current-frame arrays.
  out.LinOfOsc.assign(sz,  0LL);
  out.OscAmp.assign(sz,    0.0);
  out.OscFrq.assign(sz,    0.0);
  out.OscOfLin.assign(lsz, 0LL);  // not used by T1; kept at zero

  // Copy previous-frame state into Prv* arrays.
  out.PrvLinOfOsc.assign(sz, 0LL);
  out.PrvOscAmp.assign(sz,   0.0);
  out.PrvOscFrq.assign(sz,   0.0);

  const std::size_t prv_sz = std::min(sz, prv.LinOfOsc.size());
  for (std::size_t i = 0; i < prv_sz; ++i) {
    out.PrvLinOfOsc[i] = prv.LinOfOsc[i];
    if (i < prv.OscAmp.size()) out.PrvOscAmp[i] = prv.OscAmp[i];
    if (i < prv.OscFrq.size()) out.PrvOscFrq[i] = prv.OscFrq[i];
  }

  // [§57] Ensure hold_counts is sized and zero-extended for any new slots.
  if (hold_counts.size() < sz) hold_counts.resize(sz, 0);

  if (Nlins <= 0) {
    // No partials this frame: terminate all previously ON oscillators.
    for (std::size_t osc = 1; osc < prv_sz; ++osc) {
      if (prv.LinOfOsc[osc] > 0) ++stats.t1_terminated;
    }
    return stats;
  }

  // ── 2. Sort peaks by amplitude descending (1-based indices 1..Nlins) ─────
  std::vector<std::int64_t> peak_order(static_cast<std::size_t>(Nlins));
  std::iota(peak_order.begin(), peak_order.end(), std::int64_t{1});
  std::sort(peak_order.begin(), peak_order.end(),
    [&](std::int64_t a, std::int64_t b) {
      const double aa = (static_cast<std::size_t>(a) < LinAmp.size())
                        ? LinAmp[static_cast<std::size_t>(a)] : 0.0;
      const double bb = (static_cast<std::size_t>(b) < LinAmp.size())
                        ? LinAmp[static_cast<std::size_t>(b)] : 0.0;
      return aa > bb;  // descending
    });

  // ── 3. Per-slot claimed flag ──────────────────────────────────────────────
  std::vector<bool> osc_claimed(sz, false);

  // ── 4. DFmax at a given frequency (linear interpolation DFmax1..DFmax2) ──
  const double fc2 = (tp.Fc2 > 0.0) ? tp.Fc2 : 22050.0;
  auto dfmax_at = [&](double frq_hz) -> double {
    return tp.DFmax1 + (tp.DFmax2 - tp.DFmax1) * frq_hz / fc2;
  };

  // ── 5. Peak-claiming assignment loop ─────────────────────────────────────
  for (std::int64_t p : peak_order) {
    const std::size_t pi = static_cast<std::size_t>(p);
    if (pi >= LinFrq.size() || pi >= LinAmp.size()) continue;
    const double peak_frq = LinFrq[pi];
    const double peak_amp = LinAmp[pi];
    const double dfy = dfmax_at(peak_frq);

    // Find closest ON oscillator in prv within DFmax, not yet claimed.
    std::int64_t best_slot = 0;
    double best_dist = std::numeric_limits<double>::max();

    for (std::size_t osc = 1; osc < prv_sz; ++osc) {
      if (prv.LinOfOsc[osc] <= 0) continue;   // not ON in previous frame
      if (osc_claimed[osc])       continue;   // already claimed this frame
      const double osc_frq = (osc < prv.OscFrq.size()) ? prv.OscFrq[osc] : 0.0;
      const double dist = std::abs(osc_frq - peak_frq);
      if (dist <= dfy && dist < best_dist) {
        best_dist = dist;
        best_slot = static_cast<std::int64_t>(osc);
      }
    }

    if (best_slot > 0) {
      // Re-bind existing oscillator to this peak.
      const std::size_t bs = static_cast<std::size_t>(best_slot);
      out.LinOfOsc[bs] = p;
      out.OscFrq[bs]   = peak_frq;
      // [§64] Amplitude hysteresis: A_new = h*A_prev + (1-h)*A_peak.
      // When h=0.0 (default), A_new = A_peak (bit-identical to prior behavior).
      if (t1_amp_hysteresis > 0.0) {
        const double A_prev = (bs < prv.OscAmp.size()) ? prv.OscAmp[bs] : 0.0;
        out.OscAmp[bs] = t1_amp_hysteresis * A_prev
                       + (1.0 - t1_amp_hysteresis) * peak_amp;
      } else {
        out.OscAmp[bs] = peak_amp;
      }
      osc_claimed[bs]  = true;
      hold_counts[bs]  = 0;  // [§57] reset hold on successful claim
      ++stats.t1_assigned;
    } else {
      // Allocate the first free (previously unoccupied) slot.
      std::int64_t new_slot = 0;
      for (std::size_t osc = 1; osc < sz; ++osc) {
        const std::int64_t prv_lof = (osc < prv.LinOfOsc.size())
                                     ? prv.LinOfOsc[osc] : 0LL;
        if (!osc_claimed[osc] && prv_lof <= 0) {
          new_slot = static_cast<std::int64_t>(osc);
          break;
        }
      }
      if (new_slot > 0) {
        const std::size_t ns = static_cast<std::size_t>(new_slot);
        out.LinOfOsc[ns] = p;
        out.OscFrq[ns]   = peak_frq;
        out.OscAmp[ns]   = peak_amp;
        osc_claimed[ns]  = true;
        ++stats.t1_created;
      }
      // If no slot available (bank full), the peak is silently dropped.
    }
  }

  // ── 6. Handle unclaimed ON oscillators: hold (squelch) or terminate ─────────────────
  for (std::size_t osc = 1; osc < prv_sz; ++osc) {
    if (prv.LinOfOsc[osc] <= 0) continue;  // was not ON
    if (osc_claimed[osc]) {
      // Normally claimed this frame: hold counter already reset in step 5.
      continue;
    }
    // Not claimed.  Apply squelch-hold logic [§57].
    if (squelch_hold > 0 && hold_counts[osc] < squelch_hold) {
      // Survive: sustain last freq/amp for one more frame.
      out.LinOfOsc[osc] = prv.LinOfOsc[osc];
      out.OscFrq[osc]   = prv.OscFrq[osc];
      out.OscAmp[osc]   = prv.OscAmp[osc];
      ++hold_counts[osc];
      ++stats.t1_hold_survivals;
    } else {
      // Terminate: hold expired or hold disabled.
      // out.LinOfOsc[osc] is already 0 (step 1 init).
      // PrvOscFrq carries last known frequency for synthesis fade-out.
      if (hold_counts[osc] > 0) {
        ++stats.t1_hold_terminations;  // expired after hold
        hold_counts[osc] = 0;
      }
      ++stats.t1_terminated;
    }
  }

  // ── 7. Compute Noscs = max active slot index ──────────────────────────────
  std::int64_t noscs = 0;
  for (std::int64_t osc = MaxOscs; osc >= 1; --osc) {
    const std::size_t oi = static_cast<std::size_t>(osc);
    if (oi < out.LinOfOsc.size() && out.LinOfOsc[oi] > 0) {
      noscs = osc;
      break;
    }
  }
  out.Noscs = noscs;

  return stats;
}

} // namespace parshl
