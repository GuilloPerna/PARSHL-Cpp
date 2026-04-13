// ─────────────────────────────────────────────────────────────────────
// PARSHL — Sinusoidal Analysis/Synthesis System
// Original implementation in SAIL (1985): Julius Orion Smith III
//   Stanford CCRMA
// C++20 port (V1): Guillermo Perna (2026)
// Ref: Smith & Serra, ICMC-87 / CCRMA STAN-M-43; CMJ 14(4) 1990
// ─────────────────────────────────────────────────────────────────────
// OutPartials WAV writer — §28 / SAIL Parshl-source.txt L500–510, L1818–1824.
//
// Writes per-frame oscillator amplitude and frequency trajectories as multi-channel
// float32 WAV files at the analysis frame rate (= round(Fs / hop)).
//
//   PREFIX.amp.wav  — channel k holds OscAmp[k+1] at each analysis frame
//   PREFIX.frq.wav  — channel k holds OscFrq[k+1] at each analysis frame
//
// Both files have:
//   samplerate = round(Fs / hop)   (frame rate, not audio rate)
//   channels   = MaxOscs
//   format     = SF_FORMAT_WAV | SF_FORMAT_FLOAT  (32-bit float, IEEE)
//
// SAIL origin:
//   PROCEDURE OutPartials(MaxOscs,Nframes,Amps,Frqs,Fs,Thresh) — Parshl-source.txt L500–510
//   Interactive AmpFile/FrqFile prompts — L855–895
//   Pipeline call: OutPartials(MaxOscs,Nframes,Amps,Frqs,Fs,Thresh) — L1818–1824
//
// Invariant: flag-off behaviour (no --out-partials) is bit-identical to pre-§28.
// The tracker algorithm, synthesis, and all console output are completely unchanged.
#pragma once
#include <sndfile.h>
#include <string>
#include <vector>
#include <cstdint>

namespace parshl {

/// Writes OscAmp and OscFrq trajectories to PREFIX.amp.wav / PREFIX.frq.wav.
/// Non-copyable (owns raw SNDFILE* handles).
class PartialsWriter {
public:
    PartialsWriter() = default;
    ~PartialsWriter() noexcept { close(); }

    PartialsWriter(const PartialsWriter&)            = delete;
    PartialsWriter& operator=(const PartialsWriter&) = delete;

    /// Open PREFIX.amp.wav and PREFIX.frq.wav for writing.
    /// @param prefix        output path prefix (no extension)
    /// @param MaxOscs       oscillator bank size = number of WAV channels
    /// @param frame_rate_hz frame rate = round(Fs / hop)
    /// @return false if either file could not be opened (both are closed on failure)
    bool open(const std::string& prefix, int MaxOscs, int frame_rate_hz);

    /// Write one analysis frame. OscAmp / OscFrq are 1-based vectors
    /// (index 0 unused), size must be >= MaxOscs+1 or channels past the end
    /// are zero-filled. MaxOscs is provided for bounds-checking consistency with
    /// the caller; it is clamped to the value given to open() if larger.
    void write_frame(const std::vector<double>& OscAmp,
                     const std::vector<double>& OscFrq,
                     std::int64_t MaxOscs);

    /// Close both files. Safe to call multiple times.
    void close() noexcept;

    bool is_open() const noexcept { return sf_amp_ != nullptr; }

private:
    SNDFILE*            sf_amp_ = nullptr;
    SNDFILE*            sf_frq_ = nullptr;
    int                 MaxOscs_ = 0;
    std::vector<float>  buf_;   ///< scratch interleaved buffer (MaxOscs_ floats)
};

} // namespace parshl
