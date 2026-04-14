cmake_minimum_required(VERSION 3.20)
project(parshl_v1 CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(SNDFILE REQUIRED sndfile)
find_library(FFTW3_LIB fftw3 REQUIRED)

# ── parshl_audio_peaks — full analysis + additive synthesis driver ──────────
# V1 flags: --legacy (default, always active in V1)
add_executable(parshl_audio_peaks
  src/parshl_audio_peaks.cpp
  src/parshl_peak_claim.cpp
  src/parshl_t2_tracker.cpp
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
