/*
 * Momentum-Dependent Potential (MDYI)
 *
 * U_alpha_p_md: Core momentum-dependent mean field potential
 * From original lib/U_alpha_p_md.c
 *
 * Key physics: U(r,p) = C × Σ_j f(r_ij) / (λ²p_f0² + |p_i - p_j|²)
 *
 * Uses local_pig/g_store caching for ~40x speedup on form factors
 */

#ifndef POTENTIAL_MDYI_H
#define POTENTIAL_MDYI_H

#include "constants.h"
#include "grid_tag.h"
#include "mdyi_cache.h"

/*
 * Calculate momentum-dependent potential (MDYI mode)
 *
 * From original: lib/U_alpha_p_md.c lines 13-161
 *
 * Args:
 *   grid_tag: spatial hashing structure
 *   cache: MDYI cache (local_pig/g_store)
 *   Lx, Ly, Lz: grid cell location
 *   p_index: particle index (or 0 to clear cache)
 *   px, py, pz: momentum arrays [1..N_particles]
 *   x, y, z: position arrays [1..N_particles]
 *   grid_center_x, grid_center_y, grid_center_z: grid origins
 *   n_particles: total number of particles
 *
 * Returns:
 *   Momentum-dependent potential energy (MeV)
 *   0.0 if out of bounds or p_index=0 (cache clear signal)
 *
 * Caching behavior:
 *   - p_index=0: Clear all caches, return 0.0
 *   - First call at (Lx,Ly,Lz): Build particle list, cache form factors
 *   - Subsequent calls: Reuse cached form factors, recalculate momentum part
 */
double potential_mdyi(const GridTag *grid_tag,
                      MDYICache *cache,
                      int Lx, int Ly, int Lz,
                      int p_index,
                      const double *px, const double *py, const double *pz,
                      const double *x, const double *y, const double *z,
                      double grid_center_x, double grid_center_y, double grid_center_z,
                      int n_particles);

/*
 * Get form factor for particle at grid cell
 *
 * Helper function: computes f(r_ij) from drho_dr array
 * Original: get_R_dep() called at line 111
 *
 * Args:
 *   Lx, Ly, Lz: evaluation grid cell
 *   particle_id: which particle
 *   Lx_particle, Ly_particle, Lz_particle: particle's grid cell
 *   drho_dr: pre-computed form factor array
 *
 * Returns:
 *   Form factor value, or -1 if out of range
 */
double get_form_factor_from_drho(int Lx, int Ly, int Lz,
                                   int particle_id,
                                   int Lx_particle, int Ly_particle, int Lz_particle,
                                   const double *****drho_dr);

#endif /* POTENTIAL_MDYI_H */
