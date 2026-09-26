# LHBUU Rewrite Status

## What's Been Done (Clean, Modern Infrastructure)

### ✅ Core Architecture
- Modular design with clean headers/implementation separation
- Modern C11, no global variables
- CMake build system + manual gcc compilation
- OpenMP threading support
- SIMD vectorization pragmas
- Performance profiling and adaptive cache strategy
- Neighbor list system for efficient local calculations

### ✅ Basic Physics (Simplified - NEEDS FULL IMPLEMENTATION)
- Velocity Verlet integrator
- Grid-based density calculations
- Basic Skyrme mean-field potential
- Simplified collision dynamics
- Grid cell sorting

## What NEEDS Proper Implementation (From Original Code)

### ❌ Momentum-Dependent Potentials
**Missing from original:**
- `U_alpha_p_md.c` - Full momentum-dependent potential calculation
- `grad_V_alpha_p_md.c` - Gradient for forces
- `p_avg_alpha.c` / `p_avg_alpha_gbdmedium.c` - Local momentum averaging
- Proper double-sum over neighbors with form factors

**Current status:** Placeholder only

### ❌ Full Collision Physics
**Missing from original:**
- `scatter.c` - Complete Cugnon parametrization with differential cross sections
- `inside_sigma.c` - Proper geometric + Lorentz-contracted distance check
- `relative_velocity_ratio.c` - Medium corrections to relative velocity
- `density_states_ratio.c` - Phase space density corrections
- `p_avg_alpha_gbdmedium.c` - Momentum-dependent in-medium effects
- Proper Pauli blocking with multi-ensemble phase-space counting

**Current status:** Basic isotropic scattering, simplified blocking

### ❌ Coulomb Solver
**Missing from original:**
- Your double-grid approach (fine for nuclear, coarse for Coulomb)
- Proper boundary conditions
- FFT-based solver (FFTW integration)

**Current status:** Slow direct summation O(N²)

### ❌ Additional Potentials
**Missing from original:**
- `U_isospin.c` - Isospin asymmetry energy
- `U_alpha_sur.c` - Surface Yukawa potential
- `U_coulomb.c` - Proper Coulomb energy calculation
- Momentum-dependent isospin (MDYI)

**Current status:** Basic Coulomb only

### ❌ Initialization
**Missing from original:**
- `generate_nuclei.c` - Proper Woods-Saxon profile generation
- `resample_rogue.c` - Particle resampling
- Profile reading from files
- Damping mode for profile generation (`check_damped.c`, `damp_eom.c`)

**Current status:** Simple spherical initialization

### ❌ Diagnostics
**Missing from original:**
- `flow.c` - Flow vector analysis
- `Quad_33.c` - Quadrupole tensor
- `angular_momentum_check.c` - Detailed angular momentum
- Particle tracking by ensemble/nucleon
- Escape analysis

**Current status:** Basic energy/momentum conservation only

### ❌ Form Factors
**Missing from original:**
- Proper `form-factor.c` implementation
- Pre-computed form factor tables
- Gradient calculations for each particle (`drho_dr` array)

**Current status:** Monaghan kernel only

## Architecture Improvements Over Original

### Memory Management
- **Original:** Static compile-time arrays
- **New:** Dynamic allocation with profiling

### Threading
- **Original:** Single-threaded
- **New:** OpenMP parallelization + SIMD

### Configuration
- **Original:** Recompile for every parameter change
- **New:** Runtime configuration files

### Build System
- **Original:** Manual gcc with #include
- **New:** CMake + OpenMP detection

### Cache Optimization
- **Original:** Manual spatial sorting (brilliant for Alpha)
- **New:** Adaptive strategy based on working set size

## Next Steps (In Priority Order)

1. **Momentum-Dependent Potentials** - Core physics, expensive calculation
2. **Full Collision Module** - Port all physics from your scatter.c
3. **Coulomb FFT Solver** - Replace O(N²) with FFT
4. **Proper Initialization** - Woods-Saxon profiles
5. **Complete Diagnostics** - Flow, quadrupole, etc.

## Philosophy

Your original code represents **proven physics** validated over decades. The rewrite should:

- ✅ Modernize the software engineering (done)
- ✅ Add performance optimizations for modern hardware (done)
- ❌ **PRESERVE ALL PHYSICS** (in progress - needs completion)

The goal is not to simplify or "clean up" your physics - it's to make the same physics run faster and more conveniently on modern machines.

## Compilation

**Current (basic physics):**
```bash
gcc -o lhbuu src/*.c -Iinclude -lm -O3 -fopenmp
```

**Future (with FFTW):**
```bash
gcc -o lhbuu src/*.c -Iinclude -lm -lfftw3 -O3 -fopenmp -DUSE_FFTW
```

**Performance estimate (when complete):**
- Au+Au, 1000 ensembles, momentum-dependent
- Legacy code: ~30 minutes
- Optimized rewrite: ~10-30 seconds (60-180x faster)
