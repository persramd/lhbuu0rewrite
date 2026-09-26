/*
 * p_avg: Local Average Momentum for GBD Mode
 *
 * Implements Generalized Boltzmann Dynamics (GBD) mode where particles
 * couple to a locally-averaged momentum field. This reduces fluctuations
 * and models collective behavior in dense nuclear matter.
 *
 * Algorithm:
 * - For each grid site (i,j,k), compute weighted average momentum
 * - Weights = form factor contributions from all particles within ff_range
 * - Caches particle-to-grid mappings to avoid recomputing form factors
 * - Result stored in p_avg[ix][iy][iz][0..2] for px, py, pz
 *
 * Memory optimization:
 * - Static caching of particle indices and form factors per grid cell
 * - Cleared on demand with negative grid coordinate call
 * - Reduces form factor computations by ~10x in typical runs
 */

#ifndef P_AVG_H
#define P_AVG_H

#include "types.h"
#include "constants.h"
#include "drho_dr.h"

/*
 * Local average momentum field
 *
 * Dimensions: [nx][ny][nz][3]
 * - 3-component: [0]=px_avg, [1]=py_avg, [2]=pz_avg
 */
typedef struct {
    double ****data;    /* [ix][iy][iz][component] */
    int nx, ny, nz;     /* Grid dimensions */

    /* Cached particle contributions per grid cell */
    int ***n_cached;              /* [ix][iy][iz] = number of cached particles */
    int ****particle_indices;     /* [ix][iy][iz][n] = particle index */
    double ****form_factors;      /* [ix][iy][iz][n] = form factor weight */

    int initialized;              /* Initialization flag */
} PAvg;

/*
 * Initialize p_avg structure
 *
 * Args:
 *   nx, ny, nz: grid dimensions
 *
 * Returns:
 *   Allocated PAvg structure
 */
PAvg* p_avg_init(int nx, int ny, int nz);

/*
 * Free p_avg structure and all cached memory
 */
void p_avg_free(PAvg *pavg);

/*
 * Clear cached particle mappings
 *
 * Call this when particles move significantly or at regular intervals
 * to release memory and allow recaching
 *
 * Args:
 *   pavg: p_avg structure
 */
void p_avg_clear_cache(PAvg *pavg);

/*
 * Calculate local average momentum at grid sites
 *
 * GBD algorithm:
 *   For each grid point (ix, iy, iz):
 *     1. Find all particles within ff_range
 *     2. Weight each particle's momentum by its form factor at this point
 *     3. Normalize by total form factor: p_avg = Σ(g_i * p_i) / Σ(g_i)
 *
 * Args:
 *   pavg: p_avg structure
 *   grid: grid with spatial info
 *   particles: particle array
 *   n_particles: number of particles
 *   drho_cache: pre-computed form factors (from drho_dr)
 *   center_ix, center_iy, center_iz: optional localized update
 *     - If all zero: update entire grid
 *     - If nonzero: update only region within ±4*ff_range of center
 *
 * Mode flags:
 *   center_ix < 0: clear cache and return (no computation)
 *   center_ix = 0: full grid update
 *   center_ix > 0: localized update around (center_ix, center_iy, center_iz)
 *
 * OpenMP:
 *   Parallelized over grid sites (outer loop)
 */
void p_avg_calculate(PAvg *pavg, const Grid *grid, const Particle *particles,
                     int n_particles, const DrhoDr *drho_cache,
                     int center_ix, int center_iy, int center_iz);

/*
 * Get average momentum at grid point
 *
 * Args:
 *   pavg: p_avg structure
 *   ix, iy, iz: grid indices (1-based, as per legacy)
 *   component: 0=px, 1=py, 2=pz
 *
 * Returns:
 *   Local average momentum component (GeV)
 */
static inline double p_avg_get(const PAvg *pavg, int ix, int iy, int iz, int component)
{
    return pavg->data[ix][iy][iz][component];
}

/*
 * Get average momentum vector at particle position
 *
 * Interpolates p_avg from nearest grid cell
 *
 * Args:
 *   pavg: p_avg structure
 *   grid: grid for coordinate mapping
 *   particle: particle to query
 *   px_avg, py_avg, pz_avg: output pointers
 */
void p_avg_at_particle(const PAvg *pavg, const Grid *grid,
                       const Particle *particle,
                       double *px_avg, double *py_avg, double *pz_avg);

#endif /* P_AVG_H */
