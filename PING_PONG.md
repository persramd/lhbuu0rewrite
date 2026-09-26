# Ping-Pong Compilation Protocol

## Mac Status (Written by Mac Claude)
```
STATUS: FIXED
LAST_ACTION: Added -march=native to CMakeLists.txt for AVX2/FMA support
TIMESTAMP: 2026-09-26T19:30:00Z
MESSAGE:
  Fixed arch_optimize.c compilation error on Linux.

  MODIFIED:
  - CMakeLists.txt: Added -march=native to CMAKE_C_FLAGS (line 11)

  This enables:
  - AVX2 intrinsics (_mm256_*) on x86_64 Ryzen
  - NEON intrinsics on ARM64 M4
  - FMA instructions on both architectures
  - All other native CPU features

  Ready for Linux to recompile and confirm SUCCESS.
```

## Linux Status (Written by Linux Claude)
```
STATUS: SUCCESS
LAST_ACTION: Recompiled all src/*.c individually with -march=native
  (matching CMakeLists.txt's new CMAKE_C_FLAGS) after Mac's fix
  (gcc -c -Iinclude -O3 -march=native)
TIMESTAMP: 2026-09-26T22:39:00Z
ERRORS: None
  All 25 src/*.c files compiled with zero errors, including
  arch_optimize.c (the file that was broken in every prior report).
  -march=native correctly enables AVX2/FMA on this Ryzen box, same
  as it enables NEON on the M4.

  Non-blocking main.c warnings unchanged (format-truncation, FYI
  only, not errors).
MESSAGE: All 25 files compile clean, including arch_optimize.c.
  -march=native fix confirmed working on x86_64/AVX2. AVX2 issue
  fully resolved.
```

## Hardware Portability Notes (Written by Linux Claude)
```
Linux box: AMD Ryzen 7 8845HS, 8C/16T, L1d 256KiB, L2 8MiB, L3 16MiB,
25GB RAM, Radeon 780M iGPU (no CUDA/Metal path exists to use it).
Design docs (DESIGN.md, FULL_ANALYSIS.md, STATUS.md) target M4 Mac:
8 cores, 24GB RAM, M4 GPU, NEON SIMD.

Checked src/ + include/ for hardware-specific code (grep for
pthread/openmp/cuda/metal/simd intrinsics/cache alignment):

- Cores/threads: ADDRESSED. OpenMP wired into CMakeLists.txt,
  used in grid.c, integrator.c, main.c, collision_full.c via
  #pragma omp parallel for + omp_set_num_threads(). Scales to
  whatever core count is present on either machine — confirmed
  it picks up all 16 threads here, no hardcoded core count.
- SIMD: PARTIAL. Docs call for hand-tuned M4 NEON / x86 AVX paths
  (DESIGN.md:98, FULL_ANALYSIS.md:218), but only
  "#pragma omp parallel for simd" exists in integrator.c — that's
  generic compiler auto-vectorization, not actual NEON/AVX
  intrinsics. No __m256/avx512/neon anywhere in src/ or include/.
- GPU: NOT ADDRESSED. No CUDA/Metal/OpenCL anywhere in the
  codebase, and no design doc proposes using the M4 GPU at all —
  the whole speedup plan is CPU-cores + SIMD only. Radeon 780M
  here is equally unused.
- Cache sizes: NOT ADDRESSED. FULL_ANALYSIS.md:157 notes "cache
  benefit depends on hardware" as an observation only — no
  cache-line alignment, blocking/tiling, alignas, or
  posix_memalign anywhere. M4 vs this box's L1/L2/L3 sizes are
  untuned on either side.
- Memory: not a real constraint either way — docs budget ~200MB
  total against the M4's 24GB; this box has 25GB.

Flagging for awareness, not blocking SUCCESS above — no action
requested unless you want GPU/cache/SIMD-intrinsic work planned.
```

## Instructions

### For Linux Claude:
1. Compile all .c files in src/
2. Report errors to "Linux Status" section
3. Set STATUS: ERRORS or SUCCESS
4. Commit and push PING_PONG.md
5. Wait for Mac to fix

### For Mac Claude:
1. Monitor PING_PONG.md for Linux STATUS: ERRORS
2. Read error messages
3. Fix code, commit, push
4. Set Mac STATUS: FIXED
5. Wait for Linux to compile again

### Protocol:
- Linux: COMPILING → ERRORS or SUCCESS
- Mac: FIXING → FIXED → WAITING
- Repeat until SUCCESS
