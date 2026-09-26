# Design for Full Physics Implementation

## Momentum-Dependent Potential: The Core Challenge

### Original Algorithm (U_alpha_p_md.c)

**Calculation per particle i at grid location (Lx, Ly, Lz):**

```
U_md(i) = C * Σ_j [ f(r_ij) / (λ²p_f0² + |p_i - p_j|²) ]
```

Where:
- `f(r_ij)` = form factor (function of positions only)
- Sum over all j within form factor range (±2 cells = 5×5×5 = 125 cells)
- Typical j count: ~100 particles per grid evaluation

### Caching Strategy (Brilliance of Original)

**Problem:** Recalculating form factors is expensive
- `f(r_ij)` requires distance calculation + 3D interpolation
- Done for every particle i at every grid cell
- 4000 particles × 100 neighbors = 400,000 evaluations

**Solution:** Cache by grid cell
```c
local_pig[Lx][Ly][Lz][] = {indices of particles contributing to this cell}
g_store[Lx][Ly][Lz][]   = {pre-computed form factors}
```

**First particle at grid cell (Lx,Ly,Lz):**
1. Search neighboring 5×5×5 cells
2. Find all contributing particles
3. Calculate form factors
4. **Store both lists**

**Subsequent particles at same grid cell:**
1. Retrieve cached lists
2. Recalculate only momentum-dependent part (cheap)

**Memory management:**
- Allocate dynamically as cells are visited
- Clear all at end of force calculation (p_index=0 signal)
- Rebuild next timestep (positions changed)

### Required Data Structures

```c
typedef struct {
    int **particle_indices;      /* [nx][ny][nz] -> list of particle indices */
    float **form_factors;        /* [nx][ny][nz] -> list of form factor values */
    int *list_sizes;             /* [nx*ny*nz] -> count per cell */
    int *allocated;              /* [nx*ny*nz] -> allocation flag */
} MomentumDependentCache;
```

### Integration into Force Calculation

**Velocity Verlet requires forces twice per timestep:**

```
1. Calculate forces at t → get a(t)
2. Half-step momentum update: p(t+dt/2) = p(t) + 0.5*a(t)*dt
3. Full-step position update: r(t+dt) = r(t) + v(t+dt/2)*dt
4. Rebuild cache (positions changed)
5. Calculate forces at t+dt → get a(t+dt)
6. Half-step momentum update: p(t+dt) = p(t+dt/2) + 0.5*a(t+dt)*dt
7. Clear cache
```

**Cache lifetime:** One force calculation pass

### Performance Impact

**Without caching:**
- Form factor calculations: 4000 particles × 100 neighbors × 2 (twice/step) = 800,000
- Cost: ~10 FLOPs per form factor = 8M FLOPs
- At 50 GFLOPS: ~0.2 ms

**With caching:**
- Build cache once: 4000 × 100 = 400,000 form factors
- Reuse for ~40 particles/cell average
- Effective cost: 400,000 / 40 = 10,000 evaluations
- **40x speedup on form factors**

**Plus:** Momentum loop is vectorizable with cached data!

```c
#pragma omp simd reduction(+:sum)
for (int n = 0; n < list_size; n++) {
    int j = particle_list[n];
    float f = form_factor_cache[n];
    float dp2 = (pi.px - p[j].px)² + (pi.py - p[j].py)² + (pi.pz - p[j].pz)²;
    sum += f / (lambda_sq + dp2);
}
```

SIMD on M4: 4x speedup
Threading: 8x speedup
**Combined: ~1000x faster than naive implementation**

## Collision Physics: The Second Challenge

### Original scattered across multiple files:

**scatter.c:** Main scattering logic
- Cugnon cross section parametrization (energy-dependent)
- Lorentz boost to/from CM frame
- Differential angular distribution
- Monte Carlo acceptance

**inside_sigma.c:** Geometric check
- Closest approach distance
- Lorentz contraction in CM frame
- Geometric cross section filtering

**check_pauli.c:** Phase space blocking
- Multi-ensemble averaging
- Local Fermi momentum from density
- Occupancy counting in (r,p) space

**relative_velocity_ratio.c:** Medium corrections
- Effective mass modifications
- Density-dependent velocity

**density_states_ratio.c:** Phase space corrections
- Available states in medium

**p_avg_alpha_gbdmedium.c:** Momentum-averaged potential
- Used in cross section modifications

### Required collision data per pair:

```c
typedef struct {
    float r_closest;              /* Closest approach distance */
    float r_contracted;           /* Lorentz-contracted distance */
    float sigma_free;             /* Free NN cross section */
    float sigma_medium;           /* In-medium cross section */
    float prob_scatter;           /* Monte Carlo probability */
    float cos_theta_cm;           /* CM scattering angle */
    int pauli_blocked;            /* Blocking flag */
} CollisionData;
```

### Threading Strategy

**Problem:** Collision loop is O(N²) with state modification

**Original:** Serial loop over pairs within ensembles

**Modern approach:**
```c
#pragma omp parallel
{
    // Thread-local collision buffer
    CollisionData local_collisions[MAX_LOCAL];
    int num_local = 0;

    #pragma omp for schedule(dynamic)
    for (ensemble in ensembles) {
        for (particle pairs in ensemble) {
            if (test_collision()) {
                local_collisions[num_local++] = collision_data;
            }
        }
    }

    // Critical section: apply collisions
    #pragma omp critical
    {
        for (int i = 0; i < num_local; i++) {
            apply_collision(local_collisions[i]);
        }
    }
}
```

## Coulomb Solver: The Third Challenge

### Original approach (from comments):

**"Double grid" mentioned - likely:**
- Fine grid (dx = 1.0 fm) for nuclear forces
- Coarse grid (dx = 4.0 fm) for Coulomb
- Interpolation between grids

**Why this works:**
- Nuclear forces: short range (~2 fm) → need fine resolution
- Coulomb: long range → can use coarse grid
- FFT cost: O(N³ log N) → 64x fewer cells = huge speedup

### Modern implementation:

**Option 1: Pure FFT (FFTW)**
```c
fftw_plan_dft_r2c_3d(nx, ny, nz, rho, rho_k, FFTW_MEASURE);
// Solve φ(k) = 4πα ρ(k) / k²
fftw_plan_dft_c2r_3d(nx, ny, nz, phi_k, phi, FFTW_MEASURE);
```

**Option 2: Multigrid**
- Recursive coarse-grid relaxation
- O(N) scaling
- More complex implementation

**Option 3: Particle-Mesh Ewald (PME)**
- Split: short-range (direct) + long-range (FFT)
- Used in MD codes (GROMACS)
- Best for very large systems

**Recommendation:** Start with pure FFT, optimize later if needed

## Memory Layout for Performance

### Structure-of-Arrays (SoA) - Keep from original

**Original:**
```c
float x[N], y[N], z[N];
float px[N], py[N], pz[N];
```

**Why:** SIMD vectorization
```c
#pragma omp simd
for (int i = 0; i < N; i++) {
    x[i] += vx[i] * dt;  // Vectorizes perfectly
}
```

**Alternative (AoS):**
```c
struct Particle { float x,y,z, px,py,pz; } p[N];

for (int i = 0; i < N; i++) {
    p[i].x += p[i].vx * dt;  // Strided access, worse for SIMD
}
```

**Decision:** Keep SoA for hot loops, use AoS wrapper for convenience

### Cache-line alignment

**Original:** Natural alignment from static arrays

**Modern:**
```c
float *x = aligned_alloc(64, N * sizeof(float));  // 64-byte cache lines
```

## Implementation Priority

1. **Momentum-dependent cache system** - Enables correct force calculation
2. **Full collision module** - Core physics
3. **FFT Coulomb** - Performance bottleneck
4. **Validation against original** - Ensure physics correctness

## Testing Strategy

**Unit tests:**
- Form factor calculations
- Lorentz transformations
- Cross section parametrizations

**Integration tests:**
- Energy conservation (< 0.1%)
- Momentum conservation (< 1e-6)
- Angular momentum conservation

**Physics validation:**
- Run identical initial conditions as original
- Compare final phase space distributions
- Verify collision rates, Pauli blocking rates

**Performance benchmarks:**
- p+Ca40: Should be ~1 second for 500 steps
- Au+Au: Target <30 seconds for 500 steps

This design preserves all your physics while enabling modern optimizations.
