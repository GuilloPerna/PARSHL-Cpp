// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// OutPartials WAV writer implementation — §28 / SAIL Parshl-source.txt L500–510, L1818–1824.
//
// See parshl_partials_writer.hpp for full documentation.
#include "parshl_partials_writer.hpp"
#include <iostream>

namespace parshl {

bool PartialsWriter::open(const std::string& prefix, int MaxOscs, int frame_rate_hz) {
    close(); // ensure clean state if called more than once

    SF_INFO info{};
    info.channels   = MaxOscs;
    info.samplerate = frame_rate_hz;
    info.format     = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

    const std::string amp_path = prefix + ".amp.wav";
    const std::string frq_path = prefix + ".frq.wav";

    sf_amp_ = sf_open(amp_path.c_str(), SFM_WRITE, &info);
    if (!sf_amp_) {
        std::cerr << "[OutPartials] failed to open " << amp_path
                  << ": " << sf_strerror(nullptr) << "\n";
        return false;
    }

    sf_frq_ = sf_open(frq_path.c_str(), SFM_WRITE, &info);
    if (!sf_frq_) {
        std::cerr << "[OutPartials] failed to open " << frq_path
                  << ": " << sf_strerror(nullptr) << "\n";
        sf_close(sf_amp_);
        sf_amp_ = nullptr;
        return false;
    }

    MaxOscs_ = MaxOscs;
    buf_.assign(static_cast<std::size_t>(MaxOscs), 0.0f);

    std::cout << "[OutPartials] opened " << amp_path << " and " << frq_path
              << " (channels=" << MaxOscs << " rate=" << frame_rate_hz << " Hz)\n";
    return true;
}

void PartialsWriter::write_frame(const std::vector<double>& OscAmp,
                                 const std::vector<double>& OscFrq,
                                 std::int64_t /*MaxOscs_hint*/) {
    if (!sf_amp_ || !sf_frq_ || MaxOscs_ == 0) return;

    const std::size_t n = static_cast<std::size_t>(MaxOscs_);

    // OscAmp/OscFrq are 1-based: index 0 is unused; valid data is at [1..MaxOscs].
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t idx = i + 1; // 1-based oscillator index
        buf_[i] = (idx < OscAmp.size()) ? static_cast<float>(OscAmp[idx]) : 0.0f;
    }
    sf_writef_float(sf_amp_, buf_.data(), 1);

    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t idx = i + 1;
        buf_[i] = (idx < OscFrq.size()) ? static_cast<float>(OscFrq[idx]) : 0.0f;
    }
    sf_writef_float(sf_frq_, buf_.data(), 1);
}

void PartialsWriter::close() noexcept {
    if (sf_amp_) { sf_close(sf_amp_); sf_amp_ = nullptr; }
    if (sf_frq_) { sf_close(sf_frq_); sf_frq_ = nullptr; }
    MaxOscs_ = 0;
    buf_.clear();
}

} // namespace parshl
