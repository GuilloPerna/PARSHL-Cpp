// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port and archaeological extensions (V2): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
//(SAIL-style [lo:hi] wrappers)
#pragma once
#include <vector>
#include <stdexcept>
#include <string>
#include <cstddef>

namespace sail {

// 1D array with explicit lower/upper bounds (SAIL-style).
// Access is by the historical index (can be 0-based, 1-based, etc).
template <typename T>
class Array1D {
public:
  Array1D() = default;

  Array1D(long lo, long hi)
  : lo_(lo), hi_(hi), data_(size_from_bounds(lo, hi)) {}

  void reset(long lo, long hi) {
    lo_ = lo;
    hi_ = hi;
    data_.assign(size_from_bounds(lo, hi), T{});
  }

  [[nodiscard]] long lo() const noexcept { return lo_; }
  [[nodiscard]] long hi() const noexcept { return hi_; }
  [[nodiscard]] long size() const noexcept { return (hi_ >= lo_) ? (hi_ - lo_ + 1) : 0; }

  T& operator[](long i) {
    bounds_check(i);
    return data_[static_cast<std::size_t>(i - lo_)];
  }
  const T& operator[](long i) const {
    bounds_check(i);
    return data_[static_cast<std::size_t>(i - lo_)];
  }

  [[nodiscard]] T* raw_data() noexcept { return data_.data(); }
  [[nodiscard]] const T* raw_data() const noexcept { return data_.data(); }

private:
  long lo_ = 0;
  long hi_ = -1;
  std::vector<T> data_;

  static std::size_t size_from_bounds(long lo, long hi) {
    if (hi < lo) return 0;
    return static_cast<std::size_t>(hi - lo + 1);
  }

  void bounds_check(long i) const {
    if (i < lo_ || i > hi_) {
      throw std::out_of_range("sail::Array1D subscript out of bounds: i=" +
                              std::to_string(i) + " not in [" +
                              std::to_string(lo_) + ":" + std::to_string(hi_) + "]");
    }
  }
};

} // namespace sail
