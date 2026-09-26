/*
 * Momentum-Dependent Potential Implementation
 *
 * Direct port from lib/U_alpha_p_md.c
 */

#include "potential_mdyi.h"
#include "form_factor.h"
#include "constants.h"
#include "macros.h"
#include <stdlib.h>
#include <math.h>

/* Static normalization coefficient (lines 50-52)
 * coeff = 2 × C × λ² / ρ₀ */
static const double coeff = 2.0 * C_MDYI * LAMBDA_MDYI * LAMBDA_MDYI / RHO_0;

/*
 * Get form factor from drho_dr array
 *
 * This would normally call get_R_dep() which reads from drho_dr
 * For now, placeholder that computes form factor directly
 */
double get_form_factor_from_drho(int Lx, int Ly, int Lz,
                                   int particle_id,
                                   int Lx_particle, int Ly_particle, int Lz_particle,
                                   const double *****drho_dr)
{
    /* Calculate offset from particle's grid cell */
    int dx = Lx - Lx_particle;
    int dy = Ly - Ly_particle;
    int dz = Lz - Lz_particle;

    /* Check if within form factor range */
    if (abs(dx) > FF_RANGE || abs(dy) > FF_RANGE || abs(dz) > FF_RANGE) {
        return -1.0;
    }

    /* Convert to drho_dr array indices (0-based offset) */
    int idx_x = dx + FF_RANGE;
    int idx_y = dy + FF_RANGE;
    int idx_z = dz + FF_RANGE;

    /* Return form factor value (component [3]) */
    if (drho_dr) {
        return drho_dr[particle_id][idx_x][idx_y][idx_z][3];
    }

    /* Fallback: compute directly if drho_dr not available */
    return -1.0;
}

/*
 * Momentum-dependent potential (MDYI mode)
 *
 * Direct port from lib/U_alpha_p_md.c lines 13-161
 */
double potential_mdyi(const GridTag *grid_tag,
                      MDYICache *cache,
                      int Lx, int Ly, int Lz,
                      int p_index,
                      const double *px, const double *py, const double *pz,
                      const double *x, const double *y, const double *z,
                      double grid_center_x, double grid_center_y, double grid_center_z,
                      int n_particles)
{
    const int ff_range = FF_RANGE;
    const int max_x_grid = NX;
    const int max_y_grid = NY;
    const int max_z_grid = NZ;

    /* Normalization (line 52) */
    double pf0_hbarc = P_F0_MEV / HBARC_MEV;  /* Convert to fm^-1 */
    double norm = coeff * square(pf0_hbarc) * FF_NORM;

    /* Handle cache clear signal (lines 66-76) */
    if (p_index == 0) {
        mdyi_cache_clear(cache);
        return 0.0;
    }

    /* Make p_index positive (line 54) */
    p_index = abs(p_index);

    /* Check for boundary conditions (lines 78-81) */
    if (Lx < (ff_range + 1) || Ly < (ff_range + 1) || Lz < (ff_range + 1)) {
        return 0.0;
    }
    if (Lx > (max_x_grid - ff_range) ||
        Ly > (max_y_grid - ff_range) ||
        Lz > (max_z_grid - ff_range)) {
        return 0.0;
    }

    /* Get particle's momentum (lines 85-87) */
    double px_i = px[p_index];
    double py_i = py[p_index];
    double pz_i = pz[p_index];

    double sum = 0.0;

    /* Check if we have cached data for this cell (line 88) */
    if (!mdyi_cache_has_entry(cache, Lx, Ly, Lz)) {
        /* First evaluation at this cell: build cache (lines 91-143) */

        int search_counter = 0;
        int *pig_temp = malloc((MAX_PARTICLES + 1) * sizeof(int));
        double *g_temp = malloc((MAX_PARTICLES + 1) * sizeof(double));

        /* Search range (lines 94-99) */
        int Lx_min = Lx - ff_range;
        int Ly_min = Ly - ff_range;
        int Lz_min = Lz - ff_range;
        int Lx_max = Lx + ff_range;
        int Ly_max = Ly + ff_range;
        int Lz_max = Lz + ff_range;

        /* Loop over neighboring cells (lines 103-127) */
        for (int Lx_loop = Lx_min; Lx_loop <= Lx_max; Lx_loop++) {
            for (int Ly_loop = Ly_min; Ly_loop <= Ly_max; Ly_loop++) {
                for (int Lz_loop = Lz_min; Lz_loop <= Lz_max; Lz_loop++) {

                    /* Get particles in this cell (line 108) */
                    int m = grid_tag_count(grid_tag, Lx_loop, Ly_loop, Lz_loop);

                    for (int n = 1; n <= m; n++) {
                        int j = grid_tag_particle(grid_tag, Lx_loop, Ly_loop, Lz_loop, n);

                        /* Calculate form factor (line 111) */
                        /* This would call get_R_dep() which reads from drho_dr */
                        /* For now, compute directly */
                        int Lx_j = where_am_I(x[j], grid_center_x, INV_DX);
                        int Ly_j = where_am_I(y[j], grid_center_y, INV_DX);
                        int Lz_j = where_am_I(z[j], grid_center_z, INV_DX);

                        /* Compute form factor at evaluation point */
                        int sign[4];
                        double extras[4];
                        double g_x = form_factor_1d(x[j], Lx, grid_center_x, &sign[1], &extras[1]);
                        if (g_x <= 0) continue;

                        double g_y = form_factor_1d(y[j], Ly, grid_center_y, &sign[2], &extras[2]);
                        if (g_y <= 0) continue;

                        double g_z = form_factor_1d(z[j], Lz, grid_center_z, &sign[3], &extras[3]);
                        if (g_z <= 0) continue;

                        double f = g_x * g_y * g_z;

                        /* Store particle and form factor (lines 115-117) */
                        search_counter++;
                        pig_temp[search_counter] = j;
                        g_temp[search_counter] = f;

                        /* Calculate momentum-dependent part (lines 120-126) */
                        double delta_px = px_i - px[j];
                        double delta_py = py_i - py[j];
                        double delta_pz = pz_i - pz[j];

                        double lambda_pf0 = LAMBDA_MDYI * pf0_hbarc;
                        sum += f / (square(lambda_pf0) +
                                    square(delta_px) + square(delta_py) + square(delta_pz));
                    }
                }
            }
        }

        /* Cache the particle list if any found (lines 129-143) */
        if (search_counter > 0) {
            mdyi_cache_store(cache, Lx, Ly, Lz, search_counter, pig_temp, g_temp);
        }

        free(pig_temp);
        free(g_temp);

    } else {
        /* Cached data exists: reuse form factors (lines 146-159) */

        int search_counter = mdyi_cache_get_count(cache, Lx, Ly, Lz);

        for (int n = 1; n <= search_counter; n++) {
            int j = mdyi_cache_get_particle(cache, Lx, Ly, Lz, n);
            double f = mdyi_cache_get_form_factor(cache, Lx, Ly, Lz, n);

            /* Recalculate momentum-dependent part (lines 153-159) */
            double delta_px = px_i - px[j];
            double delta_py = py_i - py[j];
            double delta_pz = pz_i - pz[j];

            double lambda_pf0 = LAMBDA_MDYI * pf0_hbarc;
            sum += f / (square(lambda_pf0) +
                        square(delta_px) + square(delta_py) + square(delta_pz));
        }
    }

    /* Return normalized potential (implied at end) */
    return norm * sum;
}
