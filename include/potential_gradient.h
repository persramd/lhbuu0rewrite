/*
 * Momentum-Dependent Potential Gradient (MDYI)
 *
 * grad_V_alpha_p_md: Momentum-space gradient of mean field potential
 * From original lib/U_alpha_p_md.c-MDYI lines 288-478
 *
 * Key physics: ∇_p U(r,p) for momentum-dependent mean field
 *
 * The gradient is computed as:
 *   ∇_p U_i = C × Σ_{L∈ff_range(i)} f_i(L) × Σ_j f_j(L) × (p_i - p_j) / D³
 *   where D = λ²p_f0² + |p_i - p_j|²
 *
 * This is the force on particle i due to the momentum-dependent interaction.
 * Uses double-caching strategy (local_pig/g_store) for ~100x speedup.
 */

#ifndef POTENTIAL_GRADIENT_H
#define POTENTIAL_GRADIENT_H

#include "constants.h"
#include "grid_tag.h"
#include "mdyi_cache.h"
#include "drho_dr.h"

/*
 * Gradient cache structure
 *
 * Separate from MDYI potential cache because gradient needs nested loops:
 * - Outer loop: cells within ff_range of particle i
 * - Inner loop: particles contributing to each cell
 *
 * Same local_pig/g_store pattern as potential calculation
 */
typedef struct {
    MDYICache *cache;  /* Reuse MDYICache structure */
    int initialized;
} GradientCache;

/*
 * Initialize gradient cache
 *
 * Args:
 *   nx, ny, nz: grid dimensions
 *
 * Returns:
 *   Allocated GradientCache structure
 */
GradientCache* gradient_cache_init(int nx, int ny, int nz);

/*
 * Free gradient cache
 */
void gradient_cache_free(GradientCache *cache);

/*
 * Calculate momentum-space gradient of MDYI potential
 *
 * From original: lib/U_alpha_p_md.c-MDYI lines 294-478
 *
 * Computes ∇_p H where H is the total Hamiltonian minus kinetic energy.
 * This gives the force in momentum space due to momentum-dependent interactions.
 *
 * Args:
 *   grid_tag: spatial hashing structure
 *   cache: gradient cache (local_pig/g_store for each cell)
 *   drho_dr: pre-computed form factors and gradients
 *   Lx, Ly, Lz: particle's grid cell location
 *   p_index: particle index (or 0 to clear cache)
 *   px, py, pz: momentum arrays [1..N_particles]
 *   x, y, z: position arrays [1..N_particles]
 *   grid_center_x, grid_center_y, grid_center_z: grid origins
 *   n_particles: total number of particles
 *   p_grad: (output) gradient components [0..2] (0-based: x,y,z)
 *
 * Caching behavior:
 *   - p_index=0: Clear all caches, return without computing
 *   - First call at (Lx,Ly,Lz): Build particle lists, cache form factors
 *   - Subsequent calls: Reuse cached form factors, recalculate momentum part
 *
 * Boundary handling:
 *   - Requires 2×ff_range margin (double loop structure)
 *   - Returns zero gradient if too close to boundary
 *
 * Performance notes:
 *   - Double loop: outer over cells, inner over particles per cell
 *   - Form factors cached per cell (position-dependent only)
 *   - Momentum terms recalculated each time (momentum changes during timestep)
 *   - SIMD opportunity: momentum difference accumulation (see implementation)
 */
void grad_V_alpha_p_md(const GridTag *grid_tag,
                       GradientCache *cache,
                       const DrhoDr *drho_dr,
                       int Lx, int Ly, int Lz,
                       int p_index,
                       const double *px, const double *py, const double *pz,
                       const double *x, const double *y, const double *z,
                       double grid_center_x, double grid_center_y, double grid_center_z,
                       int n_particles,
                       double *p_grad);

#endif /* POTENTIAL_GRADIENT_H */
