/*
 * Momentum-Dependent Potential Gradient Implementation
 *
 * Direct port from lib/U_alpha_p_md.c-MDYI lines 294-478
 *
 * Physics:
 * --------
 * The momentum-dependent mean field creates a force in momentum space.
 * For Yukawa-type interaction: U ∝ f(r) / (λ²p_f0² + |p_i - p_j|²)
 *
 * The gradient ∇_p U determines how particle momentum evolves:
 *   dp/dt = -∇_r H + ∇_p H    (relativistic transport equation)
 *
 * Algorithm:
 * ----------
 * 1. Loop over cells L within ff_range of particle i's position
 * 2. For each cell L, sum contributions from all particles j at that cell
 * 3. Each contribution: f_i(L) × f_j(L) × (p_i - p_j) / D³
 *    where D = λ²p_f0² + |p_i - p_j|²
 *
 * Caching strategy:
 * -----------------
 * For each cell L, cache the list of particles j and their form factors f_j(L).
 * This avoids recalculating form factors which only depend on positions.
 * Momentum-dependent terms must be recalculated since p changes during timestep.
 */

#include "potential_gradient.h"
#include "form_factor.h"
#include "constants.h"
#include "macros.h"
#include "arch_optimize.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Static normalization coefficient (lines 338-347)
 * coeff = 2 × C × λ² × N_ens × dx³ / ρ₀
 * Additional factor: (-2) × (form_factor_norm)² */
static double compute_norm(void)
{
    const double pf0_hbarc = P_F0_MEV / HBARC_MEV;  /* Convert to fm^-1 */
    const double ff_norm = form_factor_norm_r(N_ENSEMBLES);

    const double coeff = 2.0 * C_MDYI * LAMBDA_MDYI * LAMBDA_MDYI
                        * N_ENSEMBLES * DX * DX * DX / RHO_0;

    return coeff * square(pf0_hbarc) * (-2.0) * square(ff_norm);
}

/*
 * Initialize gradient cache
 */
GradientCache* gradient_cache_init(int nx, int ny, int nz)
{
    GradientCache *cache = malloc(sizeof(GradientCache));
    if (!cache) return NULL;

    cache->cache = mdyi_cache_init(nx, ny, nz);
    cache->initialized = 1;

    return cache;
}

/*
 * Free gradient cache
 */
void gradient_cache_free(GradientCache *cache)
{
    if (!cache) return;

    mdyi_cache_free(cache->cache);
    free(cache);
}

/*
 * Get form factor from drho_dr array
 *
 * Helper function to extract f(x,y,z) from cached gradient array
 * Component [3] stores the form factor value
 */
static inline double get_form_factor_drho(const DrhoDr *drho_dr,
                                           int Lx_eval, int Ly_eval, int Lz_eval,
                                           int particle_id,
                                           int Lx_particle, int Ly_particle, int Lz_particle)
{
    /* Calculate offset from evaluation point to particle's cell */
    int dx = Lx_eval - Lx_particle;
    int dy = Ly_eval - Ly_particle;
    int dz = Lz_eval - Lz_particle;

    /* Check if within form factor range */
    if (abs(dx) > FF_RANGE || abs(dy) > FF_RANGE || abs(dz) > FF_RANGE) {
        return -1.0;
    }

    /* Convert to drho_dr array indices (0-based offset from -FF_RANGE to +FF_RANGE) */
    int idx_x = dx + FF_RANGE;
    int idx_y = dy + FF_RANGE;
    int idx_z = dz + FF_RANGE;

    /* Return form factor value (component [3]) */
    return drho_dr_get(drho_dr, particle_id, idx_x, idx_y, idx_z, 3);
}

/*
 * Calculate momentum-space gradient of MDYI potential
 *
 * Direct port from lib/U_alpha_p_md.c-MDYI lines 294-478
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
                       double *p_grad)
{
    const int ff_range = FF_RANGE;
    const int max_x_grid = NX;
    const int max_y_grid = NY;
    const int max_z_grid = NZ;

    /* Normalization (lines 345-347) */
    static double norm = 0.0;
    if (norm == 0.0) {
        norm = compute_norm();
    }

    /* Precompute momentum-dependent denominator coefficient */
    const double pf0_hbarc = P_F0_MEV / HBARC_MEV;
    const double lambda_pf0_sq = square(LAMBDA_MDYI * pf0_hbarc);

    /* Handle cache clear signal (lines 358-369) */
    if (p_index == 0) {
        mdyi_cache_clear(cache->cache);
        return;
    }

    /* Get particle's momentum (lines 371-373) */
    double px_i = px[p_index];
    double py_i = py[p_index];
    double pz_i = pz[p_index];

    /* Check for boundary conditions - requires 2×ff_range margin (lines 376-379) */
    if (Lx < (2*ff_range + 1) || Ly < (2*ff_range + 1) || Lz < (2*ff_range + 1)) {
        p_grad[0] = p_grad[1] = p_grad[2] = 0.0;
        return;
    }
    if (Lx > (max_x_grid - 2*ff_range) ||
        Ly > (max_y_grid - 2*ff_range) ||
        Lz > (max_z_grid - 2*ff_range)) {
        p_grad[0] = p_grad[1] = p_grad[2] = 0.0;
        return;
    }

    /* Search range around particle i (lines 380-385) */
    int Lx_min = Lx - ff_range;
    int Ly_min = Ly - ff_range;
    int Lz_min = Lz - ff_range;
    int Lx_max = Lx + ff_range;
    int Ly_max = Ly + ff_range;
    int Lz_max = Lz + ff_range;

    /* Accumulators for gradient components (line 389) */
    double grad_px = 0.0;
    double grad_py = 0.0;
    double grad_pz = 0.0;

    /* Get particle i's grid location for form factor evaluation */
    int Lx_i = where_am_I(x[p_index], grid_center_x, INV_DX);
    int Ly_i = where_am_I(y[p_index], grid_center_y, INV_DX);
    int Lz_i = where_am_I(z[p_index], grid_center_z, INV_DX);

    /* Outer loop: over cells within ff_range of particle i (lines 390-472) */
    for (int Lx_loop = Lx_min; Lx_loop <= Lx_max; Lx_loop++) {
        for (int Ly_loop = Ly_min; Ly_loop <= Ly_max; Ly_loop++) {
            for (int Lz_loop = Lz_min; Lz_loop <= Lz_max; Lz_loop++) {

                /* Get form factor g1 = f_i(Lx_loop, Ly_loop, Lz_loop) (line 393) */
                double g1 = get_form_factor_drho(drho_dr,
                                                  Lx_loop, Ly_loop, Lz_loop,
                                                  p_index,
                                                  Lx_i, Ly_i, Lz_i);

                if (g1 <= 0) continue;  /* Particle i doesn't contribute to this cell */

                /* Check if we have cached data for this cell (line 397) */
                if (!mdyi_cache_has_entry(cache->cache, Lx_loop, Ly_loop, Lz_loop)) {
                    /* First evaluation at this cell: build cache (lines 400-450) */

                    int search_counter = 0;
                    int *pig_temp = malloc((MAX_PARTICLES + 1) * sizeof(int));
                    double *g2_temp = malloc((MAX_PARTICLES + 1) * sizeof(double));

                    /* Inner search range (lines 401-406) */
                    int Lx_min2 = Lx_loop - ff_range;
                    int Ly_min2 = Ly_loop - ff_range;
                    int Lz_min2 = Lz_loop - ff_range;
                    int Lx_max2 = Lx_loop + ff_range;
                    int Ly_max2 = Ly_loop + ff_range;
                    int Lz_max2 = Lz_loop + ff_range;

                    /* Loop over neighboring cells (lines 407-434) */
                    for (int Lx_loop2 = Lx_min2; Lx_loop2 <= Lx_max2; Lx_loop2++) {
                        for (int Ly_loop2 = Ly_min2; Ly_loop2 <= Ly_max2; Ly_loop2++) {
                            for (int Lz_loop2 = Lz_min2; Lz_loop2 <= Lz_max2; Lz_loop2++) {

                                /* Get particles in this cell (line 410) */
                                int m = grid_tag_count(grid_tag, Lx_loop2, Ly_loop2, Lz_loop2);

                                for (int n = 1; n <= m; n++) {
                                    int j = grid_tag_particle(grid_tag, Lx_loop2, Ly_loop2, Lz_loop2, n);

                                    /* Get particle j's grid location */
                                    int Lx_j = where_am_I(x[j], grid_center_x, INV_DX);
                                    int Ly_j = where_am_I(y[j], grid_center_y, INV_DX);
                                    int Lz_j = where_am_I(z[j], grid_center_z, INV_DX);

                                    /* Get form factor g2 = f_j(Lx_loop, Ly_loop, Lz_loop) (line 413) */
                                    double g2 = get_form_factor_drho(drho_dr,
                                                                      Lx_loop, Ly_loop, Lz_loop,
                                                                      j,
                                                                      Lx_j, Ly_j, Lz_j);

                                    if (g2 > 0) {
                                        /* Store particle and form factor (lines 418-420) */
                                        search_counter++;
                                        pig_temp[search_counter] = j;
                                        g2_temp[search_counter] = g2;

                                        /* Calculate momentum-dependent part (lines 423-433) */
                                        double delta_px = px_i - px[j];
                                        double delta_py = py_i - py[j];
                                        double delta_pz = pz_i - pz[j];

                                        double denom = lambda_pf0_sq +
                                                      square(delta_px) +
                                                      square(delta_py) +
                                                      square(delta_pz);

                                        double temp = g2 * g1 / (denom * denom);

                                        grad_px += delta_px * temp;
                                        grad_py += delta_py * temp;
                                        grad_pz += delta_pz * temp;
                                    }
                                }
                            }
                        }
                    }

                    /* Skip caching if no particles found (line 435) */
                    if (search_counter == 0) {
                        free(pig_temp);
                        free(g2_temp);
                        continue;
                    }

                    /* Cache the particle list (lines 440-450) */
                    mdyi_cache_store(cache->cache, Lx_loop, Ly_loop, Lz_loop,
                                    search_counter, pig_temp, g2_temp);

                    free(pig_temp);
                    free(g2_temp);

                } else {
                    /* Cached data exists: reuse form factors (lines 453-471) */

                    int search_counter = mdyi_cache_get_count(cache->cache,
                                                              Lx_loop, Ly_loop, Lz_loop);

                    /* SIMD optimization opportunity: vectorize this loop
                     * - Load cached form factors g2 (position-dependent, unchanged)
                     * - Compute momentum differences (changes each timestep)
                     * - Accumulate gradient components
                     * This inner loop is hot path and benefits from NEON/AVX2 */
                    for (int n = 1; n <= search_counter; n++) {
                        int j = mdyi_cache_get_particle(cache->cache,
                                                        Lx_loop, Ly_loop, Lz_loop, n);
                        double g2 = mdyi_cache_get_form_factor(cache->cache,
                                                               Lx_loop, Ly_loop, Lz_loop, n);

                        /* Recalculate momentum-dependent part (lines 460-470) */
                        double delta_px = px_i - px[j];
                        double delta_py = py_i - py[j];
                        double delta_pz = pz_i - pz[j];

                        double denom = lambda_pf0_sq +
                                      square(delta_px) +
                                      square(delta_py) +
                                      square(delta_pz);

                        double temp = g2 * g1 / (denom * denom);

                        grad_px += delta_px * temp;
                        grad_py += delta_py * temp;
                        grad_pz += delta_pz * temp;
                    }
                }
            }
        }
    }

    /* Normalize and return gradients (lines 475-477) */
    p_grad[0] = grad_px * norm;
    p_grad[1] = grad_py * norm;
    p_grad[2] = grad_pz * norm;
}
