# Pauli Blocking Implementation Report

## Physics
Pauli blocking enforces the exclusion principle by rejecting collisions into occupied 6D phase space cells (3 position + 3 momentum dimensions), with radii derived from Fermi momentum and nuclear radius. Cold nucleus filter prevents intra-nucleus collisions before thermalization.

## Implementation Details

### Files Created
- `/Users/declan/Documents/dev/osx/lhbuu/lhbuu-rewrite/include/pauli_blocking.h` (82 lines)
- `/Users/declan/Documents/dev/osx/lhbuu/lhbuu-rewrite/src/pauli_blocking.c` (358 lines)

### Key Algorithm Components

1. **Initialization** (`pauli_blocking_init`):
   - Computes phase space radii: r_pauli, p_pauli from Fermi momentum pf0
   - Calculates max_occupancy from degeneracy factors (g=4 for spin+isospin)
   - Applies isospin-dependent corrections for neutrons/protons if D_pot != 0
   - Allocates cold nucleus interaction tracker array

2. **Cold Nucleus Filter** (lines 184-194):
   - Blocks collisions between particles from same nucleus (projectile or target)
   - Only applies before either particle has interacted with opposite nucleus
   - Prevents artificial heating of cold nuclei in early timesteps

3. **Phase Space Occupancy Check** (lines 196-353):
   - Samples n_pauli_ensembles (or all ensembles if equal to N_ensembles)
   - Uses grid acceleration to find nearby particles within r_pauli
   - For each candidate particle:
     - Fast rejection: check momentum/position component-wise before computing distance
     - Apply isospin filter (neutron vs proton) if D_pot != 0
     - Count particles in same phase space cell (within r_pauli AND p_pauli)
   - Stochastic blocking: P(block) = (count-1)/max_occupancy

4. **Isospin Dependence** (lines 95-150):
   - Separate Pauli radii for neutrons vs protons
   - Based on density asymmetry: ρ_N = (N-Z)/V_N, ρ_P = Z/V_P
   - Factors: r_neutron_iso_factor, p_neutron_iso_factor (and proton equivalents)

### Optimizations

1. **Grid Acceleration**: Uses spatial grid to limit search to cells within r_pauli
   - Computes cell bounds: Lx_low to Lx_hi (same for y, z)
   - Only checks particles in nearby cells (O(local) vs O(N))

2. **Fast Rejection**: 
   - Component-wise checks (|Δpx| > p_pauli) before computing distance²
   - Avoids expensive sqrt() calls entirely (uses squared distances)
   - Early exit when blocking_count > max_occupancy + 1

3. **Memory Efficiency**:
   - Static PauliState allocated once per simulation
   - Dynamic ensemble_selected array only during check (freed immediately)
   - Cold nucleus tracker uses int array (4 bytes/particle)

4. **OpenMP Notes**:
   - Inner loops NOT parallelized (race condition on blocking_count, early exit)
   - Collision loop in calling code should be parallelized instead
   - Thread-safe random number generator needed (TODO: replace placeholder)

### Physics Constants
- RHO0 = 0.16 fm⁻³ (nuclear saturation density)
- DEGENERACY_SPIN_ISOSPIN = 4
- DEGENERACY_SPIN = 2
- Phase factor exponent: 1/6 = 0.166667

### Exact Port from Legacy
- Algorithm matches `/Users/declan/Documents/dev/osx/lhbuu/lib/check_pauli.c` line-by-line
- Cold nucleus logic: lines 163-174 (legacy) → lines 184-194 (new)
- Isospin factors: lines 93-138 (legacy) → lines 95-150 (new)
- Grid search: lines 241-293 (legacy) → lines 269-345 (new)

### Integration Points
- Added to CMakeLists.txt (line 38)
- PauliState struct exposes all necessary state
- Functions: `pauli_blocking_init`, `pauli_blocking_check`, `pauli_blocking_mark_interaction`, `pauli_blocking_free`

### Dependencies
- `grid.h` (grid_index function for cell lookup)
- `types.h` (Particle, SimState, Grid, etc.)
- Math library (powf, sqrtf, fabsf)

### TODO for Integration
- Replace `random_number()` placeholder with actual RNG from codebase
- Call `pauli_blocking_mark_interaction()` after projectile-target collisions
- Pass parameters: n_pauli_ensembles, n_full, pf0, diff_rn_rp from config
- Parallelize collision loop in collision.c (not inner Pauli loops)

## Build Status
- Source files created and added to CMakeLists.txt
- Code ready for compilation (cmake required)
- No git commit performed (as requested)
