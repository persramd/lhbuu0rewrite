/*
 * drho_dr: Density Gradient Cache Implementation
 *
 * Direct port from lib/r_p_density.c lines 107-159
 */

#include "drho_dr.h"
#include "form_factor.h"
#include "constants.h"
#include "macros.h"
#include <stdlib.h>
#include <string.h>

/*
 * Allocate 5D array: [particle][dx][dy][dz][component]
 */
DrhoDr* drho_dr_init(int n_particles, int ff_range)
{
    DrhoDr *cache = malloc(sizeof(DrhoDr));
    cache->n_particles = n_particles;
    cache->range = ff_range;
    cache->size = 2 * ff_range + 1;

    /* Allocate 5D array */
    cache->data = malloc((n_particles + 1) * sizeof(double****));

    for (int p = 0; p <= n_particles; p++) {
        cache->data[p] = malloc(cache->size * sizeof(double***));

        for (int dx = 0; dx < cache->size; dx++) {
            cache->data[p][dx] = malloc(cache->size * sizeof(double**));

            for (int dy = 0; dy < cache->size; dy++) {
                cache->data[p][dx][dy] = malloc(cache->size * sizeof(double*));

                for (int dz = 0; dz < cache->size; dz++) {
                    cache->data[p][dx][dy][dz] = calloc(4, sizeof(double));
                }
            }
        }
    }

    return cache;
}

/*
 * Free drho_dr cache
 */
void drho_dr_free(DrhoDr *cache)
{
    if (!cache) return;

    for (int p = 0; p <= cache->n_particles; p++) {
        for (int dx = 0; dx < cache->size; dx++) {
            for (int dy = 0; dy < cache->size; dy++) {
                for (int dz = 0; dz < cache->size; dz++) {
                    free(cache->data[p][dx][dy][dz]);
                }
                free(cache->data[p][dx][dy]);
            }
            free(cache->data[p][dx]);
        }
        free(cache->data[p]);
    }
    free(cache->data);
    free(cache);
}

/*
 * Clear drho_dr array
 */
void drho_dr_clear(DrhoDr *cache)
{
    for (int p = 1; p <= cache->n_particles; p++) {
        for (int dx = 0; dx < cache->size; dx++) {
            for (int dy = 0; dy < cache->size; dy++) {
                for (int dz = 0; dz < cache->size; dz++) {
                    cache->data[p][dx][dy][dz][0] = 0.0;
                    cache->data[p][dx][dy][dz][1] = 0.0;
                    cache->data[p][dx][dy][dz][2] = 0.0;
                    cache->data[p][dx][dy][dz][3] = 0.0;
                }
            }
        }
    }
}

/*
 * Build drho_dr array for all particles
 *
 * Direct port from lib/r_p_density.c lines 107-159
 */
void drho_dr_build(DrhoDr *cache,
                   const double *x, const double *y, const double *z,
                   double grid_center_x, double grid_center_y, double grid_center_z)
{
    const int ff_range = cache->range;

    /* Loop over all particles (1-based indexing) */
    for (int k = 1; k <= cache->n_particles; k++) {
        double x_k = x[k];
        double y_k = y[k];
        double z_k = z[k];

        /* Find particle's grid cell using where_am_I macro */
        int Lr_x = where_am_I(x_k, grid_center_x, INV_DX);
        int Lr_y = where_am_I(y_k, grid_center_y, INV_DX);
        int Lr_z = where_am_I(z_k, grid_center_z, INV_DX);

        /* Grid bounds for form factor support */
        int lx_min = Lr_x - ff_range;
        int lx_max = Lr_x + ff_range;
        int ly_min = Lr_y - ff_range;
        int ly_max = Lr_y + ff_range;
        int lz_min = Lr_z - ff_range;
        int lz_max = Lr_z + ff_range;

        /* Evaluate form factor and gradient at all neighboring cells
         * Direct port from lines 142-159 */
        for (int Lx = lx_min; Lx <= lx_max; Lx++) {
            int sign_x;
            double extras_x, g_x;

            g_x = form_factor_1d(x_k, Lx, grid_center_x, &sign_x, &extras_x);

            if (g_x > 0) {
                for (int Ly = ly_min; Ly <= ly_max; Ly++) {
                    int sign_y;
                    double extras_y, g_y;

                    g_y = form_factor_1d(y_k, Ly, grid_center_y, &sign_y, &extras_y);

                    if (g_y > 0) {
                        for (int Lz = lz_min; Lz <= lz_max; Lz++) {
                            int sign_z;
                            double extras_z, g_z;

                            g_z = form_factor_1d(z_k, Lz, grid_center_z, &sign_z, &extras_z);

                            if (g_z > 0) {
                                /* Form factor value */
                                double f = g_x * g_y * g_z;

                                /* Gradient components (chain rule)
                                 * ∂f/∂x = (∂f_x/∂x) × f_y × f_z
                                 * Using extras and sign from 1D evaluation */
                                double grad_x = sgn(sign_x) * extras_x * g_y * g_z;
                                double grad_y = sgn(sign_y) * extras_y * g_z * g_x;
                                double grad_z = sgn(sign_z) * extras_z * g_x * g_y;

                                /* Store in drho_dr array
                                 * Indices are offsets from particle's cell (0-based) */
                                int idx_x = Lx - lx_min;
                                int idx_y = Ly - ly_min;
                                int idx_z = Lz - lz_min;

                                cache->data[k][idx_x][idx_y][idx_z][0] = grad_x;
                                cache->data[k][idx_x][idx_y][idx_z][1] = grad_y;
                                cache->data[k][idx_x][idx_y][idx_z][2] = grad_z;
                                cache->data[k][idx_x][idx_y][idx_z][3] = f;
                            }
                        }
                    }
                }
            }
        }
    }
}
