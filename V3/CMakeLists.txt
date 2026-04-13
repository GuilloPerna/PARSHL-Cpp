# PARSHL — Original SAIL version by J.O. Smith III (1985); C++20/23 migration by Guillermo Perna (2026).
# V3 — Experimental phase: E1–E6 improvements over V2 baseline.
cmake_minimum_required(VERSION 3.20)
project(parshl_v3 CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(SNDFILE REQUIRED sndfile)
find_library(FFTW3_LIB fftw3 REQUIRED)

# ── parshl_audio_peaks — full analysis + additive synthesis driver ──────────
# V3 features (all V2 flags, plus):
#   --oscphs-double     (E4) fractional phase preservation — RECOMMENDED
#   --closest-osc-free  (E2) closest dead-slot dispatch — RECOMMENDED
#   --sigscl-analytic   (E1) analytic SigScl = 1/sum(window)
#   --dfmax-vscale=4.0  (D7) recommended DFmax scaling
#   --dbspec-refine-amp (E3) DBspec amplitude refinement — experimental, negative Phase-A
#   --dbspec-match      (E5) DBspec candidate override — experimental, negative Phase-A
#   --grs-groups=N      (E6) GRS spectral group reduction — experimental, negative Phase-A
add_executable(parshl_audio_peaks
  src/parshl_audio_peaks.cpp
  src/parshl_partials_writer.cpp
  src/parshl_synthesize.cpp
  src/io_sndfile.cpp
  src/parshl_stft.cpp
  src/parshl_globals.cpp
  src/parshl_qinterp.cpp
  src/parshl_findpeaks.cpp
  src/parshl_dbffa.cpp
  src/fftw_wrap.cpp
  src/parshl_getwin.cpp
  src/parshl_complex_peak.cpp
  src/parshl_findpartials.cpp
  src/parshl_tracker.cpp
  src/parshl_qisort.cpp
)
target_include_directories(parshl_audio_peaks PRIVATE src ${SNDFILE_INCLUDE_DIRS})
target_link_libraries(parshl_audio_peaks PRIVATE ${SNDFILE_LIBRARIES} ${FFTW3_LIB})

# ── parshl — minimal peak-analysis CLI (no synthesis output) ────────────────
add_executable(parshl
  src/main.cpp
  src/parshl_globals.cpp
  src/parshl_qinterp.cpp
  src/parshl_findpeaks.cpp
  src/parshl_dbffa.cpp
  src/fftw_wrap.cpp
  src/parshl_complex_peak.cpp
  src/parshl_findpartials.cpp
  src/parshl_tracker.cpp
  src/parshl_qisort.cpp
)
target_include_directories(parshl PRIVATE src ${SNDFILE_INCLUDE_DIRS})
target_link_libraries(parshl PRIVATE ${SNDFILE_LIBRARIES} ${FFTW3_LIB})
