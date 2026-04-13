// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "parshl_synthesize.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace parshl {

// SAIL source: PROCEDURE Synthesize (additive synthesis via sine-table OLA)
// Parshl-source.txt L1600–L1670 (approximate range)
//
// NOTE — SAIL contained an unused fixed-point synthesis path (Parshl-source.txt L1569):
//
//   DEFINE UseFixedPoint = "FALSE";  # Set true for fast execution;
//
//   IFC UseFixedPoint THENC
//     INTEGER ARRAY SynBuf[1:Nhop];   // integer accumulators
//     INTEGER Y1, Y2, Y, Amp, Damp, Inc, Dinc, Phs;
//     DEFINE Nbits = "14";           // 14-bit fixed-point
//     DEFINE FpFp(x,y) = {((x)*(y) ASH -Nbits)};
//   ELSEC
//     REAL ARRAY SynBuf[1:Nhop];     // ← this branch is the one executed
//     REAL Y1, Y2, Y, Amp, Damp, Inc, Dinc, Phs;
//   ENDC
//
// The integer branch was never executed in the original implementation.
// PARSHL V2 implements only the floating-point branch (SAIL-faithful).
// See docs/FIXED_POINT_SYNTHESIS.md for full archaeological context.
//
// For each active oscillator, linearly interpolates amplitude and frequency
// from the previous hop to the current hop and accumulates a sine-table
// wavetable output into SynBuf[1:Nhop].  SynBuf is then overlap-added into
// the output buffer OutBuf.
//
// Key SAIL constants / relationships:
//   SinSiz  — sine table length (power of 2; SAIL OWN REAL ARRAY SinBuf)
//   Mag     — SinSiz/Fs   (converts Hz to table-samples-per-audio-sample)
//   Damp    — (OscAmp - PrvOscAmp) / Nhop   (amplitude ramp per sample)
//   Dinc    — Mag*(OscFrq - PrvOscFrq) / Nhop (frequency ramp per sample)
//   OscPhs  — INTEGER ARRAY in SAIL (L1046/L1667): assignment truncates toward
//             zero, discarding the fractional phase at each hop boundary.
//             Reproduced as cast to int (FIX C-2).
//
// Overlap-add:  SAIL L1657: OutBuf[Bp+Nsamp-1] += SynBuf[Nsamp]
void parshl_synthesize_additive(
  ParshlSynthState& ss,
  int Nhop,
  int Bp,
  int Fs,
  int Frame1,
  int MaxOscs,
  const std::vector<double>& PrvOscAmp,
  const std::vector<double>& PrvOscFrq,
  const std::vector<double>& OscAmp,
  const std::vector<double>& OscFrq,
  int NskipActive,
  std::vector<float>& OutBuf
) {
  if (Nhop <= 0 || Bp < 0 || Fs <= 0 || MaxOscs <= 0) return;

  if (Frame1 <= 1 || !ss.sine_inited) {
    // SAIL: one-time sine table initialisation (Frame1 <= 1 guard).
    // SAIL Mag <- SinSiz/Fs  (converts oscillator frequency in Hz to
    //   table-index increments per audio sample).
    const double two_pi = 2.0 * std::numbers::pi;
    const double dang = two_pi / static_cast<double>(ParshlSynthState::SinSiz);
    for (int i = 0; i < ParshlSynthState::SinSiz; ++i) {
      ss.SinBuf[static_cast<std::size_t>(i)] = static_cast<float>(std::sin(static_cast<double>(i) * dang));
    }
    ss.Mag = static_cast<double>(ParshlSynthState::SinSiz) / static_cast<double>(Fs);
    ss.sine_inited = true;
  }

  const std::size_t need_phs = static_cast<std::size_t>(MaxOscs + 1);
  if (ss.OscPhs.size() < need_phs) {
    ss.OscPhs.resize(need_phs, 0.0);
  }

  const std::size_t need_out = static_cast<std::size_t>(Bp + Nhop);
  if (OutBuf.size() < need_out) {
    OutBuf.resize(need_out, 0.0f);
  }

  std::vector<float> SynBuf(static_cast<std::size_t>(Nhop + 1), 0.0f); // 1-based like SAIL
  constexpr int sine_mask = ParshlSynthState::SinSiz - 1;

  auto synth_one_osc = [&](int Zosc) {
    if (static_cast<std::size_t>(Zosc) >= PrvOscAmp.size() ||
        static_cast<std::size_t>(Zosc) >= PrvOscFrq.size() ||
        static_cast<std::size_t>(Zosc) >= OscAmp.size() ||
        static_cast<std::size_t>(Zosc) >= OscFrq.size() ||
        static_cast<std::size_t>(Zosc) >= ss.OscPhs.size()) {
      return;
    }

    // SAIL: Damp <- (OscAmp - PrvOscAmp) / Nhop   (linear amplitude ramp)
    double Amp = PrvOscAmp[static_cast<std::size_t>(Zosc)];
    const double Damp = (OscAmp[static_cast<std::size_t>(Zosc)] - PrvOscAmp[static_cast<std::size_t>(Zosc)]) / static_cast<double>(Nhop);
    if (Amp == 0.0 && Damp == 0.0) return;

    const std::size_t Zosc_sz = static_cast<std::size_t>(Zosc);

    // SAIL: Dinc <- Mag*(OscFrq - PrvOscFrq)/Nhop  (linear frequency ramp in table steps)
    double Inc        = ss.Mag * PrvOscFrq[Zosc_sz];
    const double Dinc = ss.Mag * (OscFrq[Zosc_sz] - PrvOscFrq[Zosc_sz]) / static_cast<double>(Nhop);
    double Phs = ss.OscPhs[Zosc_sz];
    for (int Nsamp = 1; Nsamp <= Nhop; ++Nsamp) {
      // SAIL: INTEGER PROCEDURE Floor(REAL x); RETURN(x); => trunc toward zero.
      const int IntIdx    = static_cast<int>(Phs);
      const double Remain = Phs - static_cast<double>(IntIdx);
      const float Y1 = ss.SinBuf[static_cast<std::size_t>( IntIdx      & sine_mask)];
      const float Y2 = ss.SinBuf[static_cast<std::size_t>((IntIdx + 1) & sine_mask)];
      const float Y  = static_cast<float>(Y1 + Remain * (Y2 - Y1));
      SynBuf[static_cast<std::size_t>(Nsamp)] += static_cast<float>(Amp * static_cast<double>(Y));
      Amp += Damp;
      Inc += Dinc;
      Phs += Inc;
    }
    // V1 correction: fmod for full wrap robustness.
    Phs = std::fmod(Phs, static_cast<double>(ParshlSynthState::SinSiz));
    if (Phs < 0.0) Phs += static_cast<double>(ParshlSynthState::SinSiz);
    // SAIL: OscPhs is INTEGER ARRAY — truncate fractional part (audit C-2).
    ss.OscPhs[Zosc_sz] = static_cast<double>(static_cast<int>(Phs));
  };

  if (NskipActive <= 0) {
    for (int Zosc = 1; Zosc <= MaxOscs; ++Zosc) {
      synth_one_osc(Zosc);
    }
  } else {
    std::vector<int> active;
    active.reserve(static_cast<std::size_t>(MaxOscs));

    for (int Zosc = 1; Zosc <= MaxOscs; ++Zosc) {
      if (static_cast<std::size_t>(Zosc) >= PrvOscAmp.size() ||
          static_cast<std::size_t>(Zosc) >= PrvOscFrq.size() ||
          static_cast<std::size_t>(Zosc) >= OscAmp.size() ||
          static_cast<std::size_t>(Zosc) >= OscFrq.size()) {
        continue;
      }
      const double Amp = PrvOscAmp[static_cast<std::size_t>(Zosc)];
      const double Damp = (OscAmp[static_cast<std::size_t>(Zosc)] - PrvOscAmp[static_cast<std::size_t>(Zosc)]) / static_cast<double>(Nhop);
      if (Amp == 0.0 && Damp == 0.0) continue;
      active.push_back(Zosc);
    }

    std::sort(active.begin(), active.end(), [&](int a, int b) {
      const double fa = PrvOscFrq[static_cast<std::size_t>(a)];
      const double fb = PrvOscFrq[static_cast<std::size_t>(b)];
      if (fa < fb) return true;
      if (fa > fb) return false;
      return a < b;
    });

    const std::size_t skip = static_cast<std::size_t>(std::min<int>(NskipActive, static_cast<int>(active.size())));
    for (std::size_t k = skip; k < active.size(); ++k) {
      synth_one_osc(active[k]);
    }
  }

  // SAIL L1657: FOR Nsamp <- 1 STEP 1 UNTIL Nhop DO OutBuf[Bp+Nsamp-1] <- OutBuf[Bp+Nsamp-1] + SynBuf[Nsamp]
  // Overlap-add accumulation into output buffer.
  for (int Nsamp = 1; Nsamp <= Nhop; ++Nsamp) {
    OutBuf[static_cast<std::size_t>(Bp + Nsamp - 1)] += SynBuf[static_cast<std::size_t>(Nsamp)];
  }
}

} // namespace parshl
