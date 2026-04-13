// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "fftw_wrap.hpp"

namespace parshl {

FftwR2C::FftwR2C(int n) : N(n), in(n), out(n/2 + 1) {
  plan = fftw_plan_dft_r2c_1d(
    N,
    in.data(),
    out.data(),
    FFTW_ESTIMATE
  );
  if (!plan) throw std::runtime_error("FFTW: failed to create plan");
}

FftwR2C::~FftwR2C() {
  if (plan) fftw_destroy_plan(plan);
}

void FftwR2C::exec() {
  fftw_execute(plan);
}

FftwC2R::FftwC2R(int n) : N(n), in(n/2 + 1), out(n) {
  plan = fftw_plan_dft_c2r_1d(
    N,
    in.data(),
    out.data(),
    FFTW_ESTIMATE
  );
  if (!plan) throw std::runtime_error("FFTW: failed to create inverse plan");
}

FftwC2R::~FftwC2R() {
  if (plan) fftw_destroy_plan(plan);
}

void FftwC2R::exec() {
  fftw_execute(plan);
}

} // namespace parshl
