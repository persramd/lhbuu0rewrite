# LHBUU Rewrite - Current Status
**Date:** 2026-09-28  
**Repository:** https://github.com/persramd/lhbuu0rewrite.git  
**Latest Commit:** 595e912

## Build Status
- **M4 Mac:** ✓ Compiles (not yet linked)
- **AMD Ryzen Linux:** ⏳ Testing full link (duplicate symbol fix)

## Architecture Support
| Platform | CPU | SIMD | GPU | Status |
|----------|-----|------|-----|--------|
| M4 Mac | ARM64 | NEON | Metal | Infrastructure ready |
| AMD Ryzen | x86_64 | AVX2 | HIP (planned) | CPU only |

## Physics Modules Implemented

### Core Infrastructure ✓
- Form factors (Monaghan M3 cubic spline)
- Grid management (64³ cells, spatial hashing)
- MDYI potential (momentum-dependent Yukawa)
- Cache systems (drho_dr, mdyi_cache, grid_tag)

### Collision Physics ✓
- Cross sections (in-medium corrections)
- Lorentz transformations
- Pauli blocking (6D phase space + cold nucleus filter)
- Inside sigma calculation

### Forces & Integration ✓
- Half-impulse forces (nuclear, Coulomb stubs)
- Potential gradient (∇_p V for MDYI)
- Modified Velocity Verlet integrator
- GBD mode (p_avg_alpha)

### Initialization ✓
- Woods-Saxon nuclear profiles
- Medium corrections (density of states, velocity ratio)

## Pending Work

### Physics Completeness
- [ ] Charge density grid (for isospin forces)
- [ ] Surface Yukawa potential
- [ ] Wire GBD mode into integrator
- [ ] Complete MDYI grid updates

### Main Loop
- [ ] Initialize simulation state
- [ ] Time stepping with all forces
- [ ] Diagnostics output
- [ ] Validation against 1997 code

### GPU Acceleration (On Hold)
- [x] Metal infrastructure (M4)
- [x] Grid density Metal kernel
- [ ] HIP backend (Ryzen)
- [ ] Additional kernels (forces, FFT)
- [ ] Performance benchmarks

## Recent Fixes (Ping-Pong Loop)
1. Missing stdint.h → Added
2. AVX2 intrinsics fail → Added -march=native
3. OBJC breaks Linux cmake → Conditional enable_language(OBJC)
4. Duplicate collision_process_all → Removed collision.c

## Performance Target
**Goal:** 100-1000x speedup vs 1997 DEC Alpha code  
**Strategy:** OpenMP + SIMD (CPU), Metal/HIP (GPU), architecture-specific optimization

## Token Usage
~104k/200k (52% of session limit)
