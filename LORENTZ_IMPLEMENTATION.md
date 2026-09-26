# LHBUU Lorentz Transformation Implementation

## Overview

Implemented complete Lorentz transformation utilities for collision physics in the LHBUU rewrite. These functions provide boost transformations to/from the center-of-mass (CM) frame, which are essential for computing nucleon-nucleon collision cross sections and scattering dynamics.

## Source Origins

### Original Code Locations

1. **lorentz_boostnntoNN.c**
   - Path: `/Users/declan/Documents/dev/osx/lhbuu/lib/lorentz_boostnntoNN.c`
   - Purpose: High-level boost function for ensemble transformation
   - Key contribution: Overall boost formula structure

2. **collisions.c**
   - Path: `/Users/declan/Documents/dev/osx/lhbuu/lib/collisions.c`
   - Functions ported:
     - `lorentz()`: Calculate CM parameters from particle pair
     - `root_s()`: Invariant mass computation
     - Boost formula implementation from lines 220-365
   - Key contribution: Detailed boost formulas and parameter calculations

## Implementation Files

### Header: `/Users/declan/Documents/dev/osx/lhbuu/lhbuu-rewrite/include/lorentz.h`

Defines the public API with:

- **LorentzParams structure**: Contains CM frame parameters
  - `cm_energy`: Invariant mass (sqrt(s))
  - `beta_x, beta_y, beta_z`: Velocity components of CM frame
  - `beta`: Magnitude of velocity |beta|
  - `gamma`: Lorentz factor γ = 1/√(1-β²)
  - `gamma_m1_over_beta2`: Pre-computed (γ-1)/β² for numerical stability

- **Core Functions**:
  - `lorentz_params_from_pair()`: Calculate CM parameters from two particles
  - `lorentz_boost_to_cm()`: Transform particle momentum to CM frame
  - `lorentz_boost_from_cm()`: Inverse transformation back to lab frame

- **Utility Functions**:
  - `lorentz_gamma()`: Compute Lorentz factor
  - `lorentz_velocity()`: Calculate velocity from momentum
  - `lorentz_energy()`: Compute total relativistic energy
  - `lorentz_invariant_mass()`: Calculate sqrt(s) for two-particle system

### Implementation: `/Users/declan/Documents/dev/osx/lhbuu/lhbuu-rewrite/src/lorentz.c`

Contains the complete implementation using double-precision internal arithmetic for numerical stability while interfacing with float-based Particle structures.

## Key Formulas

### Center-of-Mass Velocity

The velocity of the CM frame is computed from total momentum and energy:

```
β_i = p_total_i / E_total
E_total = E1 + E2
p_total = p1 + p2
```

### Lorentz Factor

```
γ = 1 / √(1 - β²)
```

### Invariant Mass (CM Energy)

```
s = (E1 + E2)² - |p1 + p2|²
√s = CM energy
```

### Boost Transformation

The standard relativistic boost formula transforms a particle's 4-momentum when moving to a frame with velocity -β:

```
p'_coefficient = (γ-1)/β² * (p·β) - γ*E
p' = p + p'_coefficient * β
```

For the inverse transformation (CM → lab), the sign of the energy term flips:

```
p'_coefficient = (γ-1)/β² * (p·β) + γ*E
p' = p + p'_coefficient * β
```

**Key advantage**: Using (γ-1)/β² is numerically stable even for very small β, avoiding division-by-zero issues.

## Numerical Characteristics

### Units

- Momentum: GeV/c (from Particle structure)
- Energy: GeV (relativistic)
- Velocity: Units of c
- Mass: GeV/c² (NUCLEON_MASS = 0.938919 GeV from types.h)

### Precision

- Internal calculations: Double precision (C `double`)
- Particle data: Single precision (C `float`)
- Conversion happens at function boundaries

### Stability Features

1. **Pre-computed (γ-1)/β² factor**: Avoids numerical issues when β is small
2. **Direct use of sqrt()**: Avoids power function (faster, more stable)
3. **Sanity checks**: Validates β < 1 before computing transformations
4. **Error handling**: Returns -1 on invalid input (NULL pointers, unphysical parameters)

## Testing

### Test Results

Two comprehensive test programs validate the implementation:

#### Test 1: Static CM (Zero Velocity CM Frame)
- Particles with equal and opposite momentum
- Beta = 0 (particles already in CM frame)
- Results:
  - CM energy computed correctly
  - Round-trip boost accuracy: 1e-6 level
  - Status: PASS

#### Test 2: Dynamic CM (Non-zero Velocity CM Frame)
- Realistic collision scenario: projectile + target
- CM velocity: β = 0.416c (γ = 1.0996)
- Results:
  - CM energy: 1.986 GeV (physically reasonable)
  - Total momentum in CM: < 1e-8 GeV (excellent)
  - Round-trip errors: < 1e-6 GeV
  - Status: PASS

### Compilation

```bash
# Compile object file
gcc -c -I./include -o lorentz.o src/lorentz.c -Wall -Wextra -O3 -lm

# Link with test program
gcc test_lorentz.c lorentz.o -o test_lorentz -lm

# Run
./test_lorentz
```

**Output**: Successfully compiles to Mach-O 64-bit object (arm64 architecture)

## Integration with Build System

### CMakeLists.txt Update

Added `src/lorentz.c` to SOURCES list:

```cmake
set(SOURCES
    ...
    src/collision.c
    src/lorentz.c
    ...
)
```

The lorentz module is automatically compiled as part of the main lhbuu executable.

## Dependencies

### Header Files
- `types.h`: Particle structure definition
- `macros.h`: square(), mag3(), mag4() macros
- `<math.h>`: Standard math library (sqrt, etc.)
- `<stddef.h>`: NULL definition

### External Constants
- `NUCLEON_MASS` (from types.h): 0.938919 GeV

## Physical Interpretation

### Boost Direction

The boost functions implement transformation to/from the **two-particle rest frame** (binary CM frame), which is essential for collision calculations:

1. **To CM frame**: Transform initial-state particles
   - Particles appear head-on in CM frame
   - Enables simple scattering angle definitions

2. **From CM frame**: Transform final-state particles
   - Scattered particles returned to lab frame
   - Pauli blocking checks performed in final momenta

### Numerical Validation

The near-perfect round-trip accuracy (< 1e-6 relative error) confirms:
- Correct boost formula implementation
- Proper handling of gamma factor
- Accurate energy-momentum calculations
- Suitability for collision dynamics where momentum conservation must be maintained to high precision

## Future Enhancements

Potential extensions:
1. **Many-body boosts**: Transform entire ensemble to CM frame
2. **Spatial transformations**: Include coordinate boosting for Lorentz contraction
3. **Diagnostic output**: Optional printing of boost parameters for debugging
4. **Batch operations**: Vectorized boost for multiple particles
5. **Special cases**: Fast paths for collinear collisions (β_x = 0, β_y = 0)

## References

- Original code: `/Users/declan/Documents/dev/osx/lhbuu/lib/`
- Design document: `/Users/declan/Documents/dev/osx/lhbuu/lhbuu-rewrite/DESIGN.md`
- Classical special relativity: Goldstein, Poole, Safko; Peskin & Schroeder
