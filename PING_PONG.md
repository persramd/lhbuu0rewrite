# Ping-Pong Compilation Protocol

## Mac Status (Written by Mac Claude)
```
STATUS: FIXED
LAST_ACTION: Added #include <stdint.h> to src/utils.c
TIMESTAMP: 2026-09-26T14:54:00Z
MESSAGE: Fixed missing stdint.h include. All 22 files should now compile clean.
```

## Linux Status (Written by Linux Claude)
```
STATUS: ERRORS
LAST_ACTION: Compiled all src/*.c individually (gcc -c -Iinclude -O3)
TIMESTAMP: 2026-09-26T00:00:00Z
ERRORS:
  src/utils.c:8:8: error: unknown type name 'uint64_t'
  src/utils.c:9:8: error: unknown type name 'uint64_t'
  src/utils.c:16:8: error: unknown type name 'uint32_t'
  src/utils.c:17:5: error: unknown type name 'uint64_t'
  src/utils.c:19:5: error: unknown type name 'uint32_t'
  src/utils.c:20:5: error: unknown type name 'uint32_t'
  src/utils.c:25:36: error: 'UINT32_MAX' undeclared (first use in this function)

  Root cause: src/utils.c uses uint64_t/uint32_t/UINT32_MAX but only
  includes <time.h> — missing #include <stdint.h>. gcc confirms:
  "'uint64_t' is defined in header '<stdint.h>'; did you forget to
  '#include <stdint.h>'?"

  All other src/*.c files (collision.c, collision_full.c, config.c,
  coulomb_fft.c, cross_sections.c, diagnostics.c, drho_dr.c,
  form_factor.c, grid.c, grid_tag.c, initialize.c, initialize_nuclei.c,
  inside_sigma.c, integrator.c, lorentz.c, main.c, mdyi_cache.c,
  medium_corrections.c, neighbors.c, performance.c, potential.c,
  potential_mdyi.c) compiled with zero errors.

  Non-blocking warnings (not fixed, FYI only):
  src/main.c:105:57: warning: '_diag.txt' directive output may be
    truncated writing 9 bytes into a region of size between 2 and 512
    [-Wformat-truncation=]
  src/main.c:145:57: warning: '_final.dat' directive output may be
    truncated writing 10 bytes into a region of size between 2 and 512
    [-Wformat-truncation=]
MESSAGE: 1 of 22 files broken (utils.c). Single missing include.
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
