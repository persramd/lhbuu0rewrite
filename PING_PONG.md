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
STATUS: ERRORS
LAST_ACTION: Pulled GPU infrastructure push (07c4b7c: GPU_INFRASTRUCTURE.md,
  METAL_GPU_IMPLEMENTATION.md, gpu_compute.h, gpu_compute_metal.m,
  gpu_compute_cpu.c, Metal shaders, tests). Compiled all src/*.c
  individually (gcc -c -Iinclude -O3 -march=native), then ran a REAL
  `cmake .` configure to test the full build (not just per-file compiles).
TIMESTAMP: 2026-09-26T23:05:00Z
ERRORS:
  Per-file compiles: all 26 src/*.c files (25 previous + new
  gpu_compute_cpu.c) compile with zero errors. gpu_compute_cpu.c is
  clean, no issues.

  BUT the actual `cmake .` configure step FAILS outright on Linux:

  CMakeLists.txt:2 now reads `project(LHBUU C OBJC)` — this
  unconditionally declares Objective-C as a project language. CMake
  tests the OBJC compiler at configure time, before any of the
  GPU_ENABLED/APPLE conditional logic (lines 13-43) even runs. This
  Linux box's gcc has no Objective-C front-end (`cc1obj` missing), so
  configure aborts immediately:

    CMake Error at .../CMakeTestOBJCCompiler.cmake:64 (message):
      The Objective-C compiler "/usr/bin/cc" is not able to compile a
      simple test program.
      cc: fatal error: cannot execute 'cc1obj': execvp: No such file
      or directory
    CMake Error: CMake will not be able to correctly generate this
    project.
    -- Configuring incomplete, errors occurred!

  Root cause: `project(...)` language declarations are checked
  up-front regardless of whether any OBJC source ends up in the
  actual SOURCES list. Even though gpu_compute_metal.m is correctly
  excluded from SOURCES on Linux (GPU_ENABLED is FALSE via the
  UNIX-AND-NOT-APPLE branch, line 31-35, which falls back to
  gpu_compute_cpu.c), the mere presence of `OBJC` in the top-level
  `project()` call is enough to break configure here.

  This is a NEW class of error the old gcc-per-file loop can't catch
  — it only shows up when actually running `cmake .`. Recommend
  Linux Claude's compile step going forward also run a real `cmake .`
  configure (not just gcc -c per file) to catch build-system-level
  issues like this.
MESSAGE: All 26 individual .c files compile clean, including the new
  gpu_compute_cpu.c CPU fallback. However the real build is broken:
  `cmake .` fails to configure on Linux because `project(LHBUU C OBJC)`
  requires an Objective-C compiler unconditionally, and this box
  doesn't have one. Needs the OBJC language declaration made
  conditional (e.g. only add OBJC via `enable_language(OBJC)` inside
  the `if(APPLE)` block, not in the top-level project() call).
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
