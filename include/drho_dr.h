/*
 * drho_dr: Density Gradient Cache
 *
 * Pre-computed form factors and gradients for each particle
 * Enables fast force calculation without repeated form factor evaluation
 *
 * Philosophy from original:
 * - Calculate once per timestep (after position update)
 * - Reuse many times during force calculation
 * - Sparse storage: only ±ff_range cells around particle
 */

#ifndef DRHO_DR_H
#define DRHO_DR_H

#include "constants.h"

/*
 * drho_dr array structure
 *
 * Dimensions: [particle][Δx][Δy][Δz][component]
 * - particle: particle index (1..N_t)
 * - Δx,Δy,Δz: offset from particle's grid cell (0..2*FF_RANGE)
 * - component: [0]=∂f/∂x, [1]=∂f/∂y, [2]=∂f/∂z, [3]=f(x,y,z)
 *
 * Memory: MAX_PARTICLES × (2*FF_RANGE+1)³ × 4 × 8 bytes
 *         20000 × 5³ × 4 × 8 = 80 MB (manageable)
 */
typedef struct {
    double *****data;  /* [particle][dx][dy][dz][component] */
    int n_particles;
    int range;         /* FF_RANGE */
    int size;          /* 2*range + 1 */
} DrhoDr;

/*
 * Initialize drho_dr cache
 *
 * Args:
 *   n_particles: total number of particles
 *   ff_range: form factor range (±cells)
 *
 * Returns:
 *   Allocated DrhoDr structure
 */
DrhoDr* drho_dr_init(int n_particles, int ff_range);

/*
 * Free drho_dr cache
 */
void drho_dr_free(DrhoDr *cache);

/*
 * Build drho_dr array for all particles
 *
 * Computes form factors and gradients for each particle at all
 * neighboring grid cells within ff_range
 *
 * Args:
 *   cache: drho_dr structure
 *   x, y, z: particle positions [1..n_particles]
 *   grid_center_x, grid_center_y, grid_center_z: grid origins
 *
 * Called once per timestep after position update
 */
void drho_dr_build(DrhoDr *cache,
                   const double *x, const double *y, const double *z,
                   double grid_center_x, double grid_center_y, double grid_center_z);

/*
 * Clear drho_dr array (zero all entries)
 *
 * Called before rebuilding
 */
void drho_dr_clear(DrhoDr *cache);

/*
 * Access drho_dr value
 *
 * Args:
 *   cache: drho_dr structure
 *   particle: particle index (1-based)
 *   dx, dy, dz: grid offset (0..2*FF_RANGE)
 *   component: 0=∂f/∂x, 1=∂f/∂y, 2=∂f/∂z, 3=f
 *
 * Returns:
 *   Cached value
 */
static inline double drho_dr_get(const DrhoDr *cache, int particle,
                                  int dx, int dy, int dz, int component)
{
    return cache->data[particle][dx][dy][dz][component];
}

/*
 * Set drho_dr value
 *
 * Used internally during build
 */
static inline void drho_dr_set(DrhoDr *cache, int particle,
                                int dx, int dy, int dz, int component,
                                double value)
{
    cache->data[particle][dx][dy][dz][component] = value;
}

#endif /* DRHO_DR_H */
