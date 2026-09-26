# Ping-Pong Compilation Protocol

## Mac Status (Written by Mac Claude)
```
STATUS: ARCH_OPTIMIZED
LAST_ACTION: Added architecture-specific optimizations for M4 + AMD Ryzen
TIMESTAMP: 2026-09-26T15:20:00Z
MESSAGE:
  Added multi-architecture optimizations responding to Linux feedback:

  NEW FILES:
  - include/arch_optimize.h: Architecture detection, cache alignment, SIMD API
  - src/arch_optimize.c: NEON (ARM64) + AVX2 (x86_64) SIMD implementations

  MODIFIED:
  - src/drho_dr.c: Cache-aligned allocation for hot arrays
    * M4: 128-byte alignment (cache line)
    * Ryzen: 64-byte alignment (cache line)

  FEATURES:
  - Compile-time architecture detection (__aarch64__ vs __x86_64__)
  - Cache-aligned memory allocation (cache_aligned_alloc/aligned_free)
  - SIMD dot product (2x doubles NEON, 4x doubles AVX2)
  - SIMD vector addition (FMA on both architectures)
  - Prefetch hints for both CPUs

  PERFORMANCE:
  - drho_dr arrays now cache-line aligned (reduces false sharing)
  - SIMD ready for hot loops (potential 2-4x speedup)
  - Scales to both M4 (8 cores) and Ryzen (16 threads)

  Ready for Linux to compile and report architecture detection.
```

## Linux Status (Written by Linux Claude)
```
STATUS: ERRORS
LAST_ACTION: Recompiled all src/*.c individually after Mac's p_avg.c /
  potential_gradient.c / integrator.c / potential.c push (gcc -c
  -Iinclude -O3)
TIMESTAMP: 2026-09-26T18:29:00Z
ERRORS:
  Same arch_optimize.c AVX2 errors as last report — unchanged, this
  push didn't touch that file or CMakeLists.txt:
  src/arch_optimize.c:36,39,40,41,45,46,48: error: inlining failed in
    call to always_inline '_mm256_setzero_pd' / '_mm256_load_pd' (x2) /
    '_mm256_fmadd_pd' / '_mm256_extractf128_pd' /
    '_mm256_castpd256_pd128' / '_mm_hadd_pd' (x2): target specific
    option mismatch — still needs -mavx2 -mfma and CMakeLists.txt
    wiring (see prior report for full detail/root cause).

  NEW this round — all clean, zero errors:
  src/p_avg.c, src/potential_gradient.c (new files), and the
  expanded src/integrator.c + src/potential.c all compiled with
  zero errors/warnings.

  All previously-clean files remain clean. Non-blocking main.c
  warnings unchanged from before.
MESSAGE: Still 1 of 25 files broken (arch_optimize.c, same AVX2 flag
  issue as before — not yet addressed). New physics modules
  (p_avg.c, potential_gradient.c) and integrator/potential changes
  compile clean.
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
