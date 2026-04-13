// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port, archaeological extensions (V2), and experimental
//   improvements (V3): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
#include "io_sndfile.hpp"
#include <sndfile.h>
#include <iostream>
#include <vector>
#include <cstdint>

namespace parshl {

bool ReadMonoAsDouble(const std::string& path, std::vector<double>& out, int& fs) {
  SF_INFO info{};
  SNDFILE* sf = sf_open(path.c_str(), SFM_READ, &info);
  if (!sf) {
    std::cerr << "libsndfile: cannot open input: " << path << " : " << sf_strerror(nullptr) << "\n";
    return false;
  }

  fs = info.samplerate;
  const int channels = info.channels;
  if (channels <= 0) {
    std::cerr << "libsndfile: invalid channel count\n";
    sf_close(sf);
    return false;
  }

  const sf_count_t frames = info.frames;
  std::vector<double> interleaved(static_cast<size_t>(frames) * static_cast<size_t>(channels));
  sf_count_t got = sf_readf_double(sf, interleaved.data(), frames);
  if (got != frames) {
    std::cerr << "libsndfile: short read: got " << got << " of " << frames << " frames\n";
  }
  sf_close(sf);

  out.resize(static_cast<size_t>(got));
  // SANDI (DEC-10) delivered samples as 16-bit integers [-32768, 32767].
  // libsndfile normalizes to [-1, 1]. Scale factor restores original range
  // so all downstream code (SigScl, LinAmp) behaves as in the SAIL original.
  // Take channel 0 (archaeological simple choice).
  for (sf_count_t i = 0; i < got; ++i) {
    out[static_cast<size_t>(i)] = interleaved[static_cast<size_t>(i) * channels + 0] * 32768.0;
  }
  return true;
}

bool WriteMonoWavFloat32(const std::string& path, const std::vector<double>& in, int fs) {
  SF_INFO info{};
  info.channels = 1;
  info.samplerate = fs;
  info.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

  SNDFILE* sf = sf_open(path.c_str(), SFM_WRITE, &info);
  if (!sf) {
    std::cerr << "libsndfile: cannot open output: " << path << " : " << sf_strerror(nullptr) << "\n";
    return false;
  }

  std::vector<float> tmp(in.size());
  for (size_t i = 0; i < in.size(); ++i) tmp[i] = static_cast<float>(in[i]);

  sf_count_t wrote = sf_writef_float(sf, tmp.data(), static_cast<sf_count_t>(tmp.size()));
  if (wrote != static_cast<sf_count_t>(tmp.size())) {
    std::cerr << "libsndfile: short write: wrote " << wrote << " of " << tmp.size() << " frames\n";
  }
  sf_close(sf);
  return true;
}

} // namespace parshl
