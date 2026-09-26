/*
 * p_avg: Local Average Momentum for GBD Mode
 *
 * Direct port from lib/p_avg_alpha.c with modernization:
 * - C11, double precision
 * - OpenMP parallelization over grid sites
 * - Cache-aligned memory for SIMD-friendly data layout
 * - Uses modern grid structure and drho_dr cache
 *
 * Physics preserved exactly:
 * - Form factor weighted momentum averaging
 * - Particle caching optimization
 * - Localized vs full grid update modes
 */

#include "p_avg.h"
#include "grid.h"
#include "constants.h"
#include "macros.h"
#include "arch_optimize.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * Helper: get form factor contribution from drho_dr cache
 *
 * This replaces the get_R_dep() call in legacy code
 *
 * Args:
 *   drho_cache: pre-computed form factors
 *   particle_idx: particle index (1-based)
 *   grid_ix, grid_iy, grid_iz: target grid cell
 *   particle_ix, particle_iy, particle_iz: particle's home grid cell
 *
 * Returns:
 *   Form factor value g (positive if particle contributes, 0 otherwise)
 */
static inline double get_form_factor(const DrhoDr *drho_cache,
                                      int particle_idx,
                                      int grid_ix, int grid_iy, int grid_iz,
                                      int particle_ix, int particle_iy, int particle_iz)
{
    /* Calculate offset in drho_dr coordinates */
    int alpha_x = (grid_ix - particle_ix) + FF_RANGE;
    int alpha_y = (grid_iy - particle_iy) + FF_RANGE;
    int alpha_z = (grid_iz - particle_iz) + FF_RANGE;

    /* Component [3] is the form factor f(r) */
    return drho_dr_get(drho_cache, particle_idx, alpha_x, alpha_y, alpha_z, 3);
}

/*
 * Allocate p_avg structure
 */
PAvg* p_avg_init(int nx, int ny, int nz)
{
    PAvg *pavg = malloc(sizeof(PAvg));
    pavg->nx = nx;
    pavg->ny = ny;
    pavg->nz = nz;
    pavg->initialized = 1;

    /* Allocate main data array [ix][iy][iz][3] */
    pavg->data = malloc((nx + 1) * sizeof(double***));
    for (int i = 0; i <= nx; i++) {
        pavg->data[i] = malloc((ny + 1) * sizeof(double**));
        for (int j = 0; j <= ny; j++) {
            pavg->data[i][j] = malloc((nz + 1) * sizeof(double*));
            for (int k = 0; k <= nz; k++) {
                pavg->data[i][j][k] = cache_aligned_alloc(3 * sizeof(double));
                memset(pavg->data[i][j][k], 0, 3 * sizeof(double));
            }
        }
    }

    /* Allocate cache arrays */
    pavg->n_cached = malloc((nx + 1) * sizeof(int**));
    pavg->particle_indices = malloc((nx + 1) * sizeof(int***));
    pavg->form_factors = malloc((nx + 1) * sizeof(double***));

    for (int i = 0; i <= nx; i++) {
        pavg->n_cached[i] = malloc((ny + 1) * sizeof(int*));
        pavg->particle_indices[i] = malloc((ny + 1) * sizeof(int**));
        pavg->form_factors[i] = malloc((ny + 1) * sizeof(double**));

        for (int j = 0; j <= ny; j++) {
            pavg->n_cached[i][j] = calloc(nz + 1, sizeof(int));
            pavg->particle_indices[i][j] = calloc(nz + 1, sizeof(int*));
            pavg->form_factors[i][j] = calloc(nz + 1, sizeof(double*));
        }
    }

    return pavg;
}

/*
 * Free p_avg structure
 */
void p_avg_free(PAvg *pavg)
{
    if (!pavg) return;

    /* Free cached data */
    for (int i = 0; i <= pavg->nx; i++) {
        for (int j = 0; j <= pavg->ny; j++) {
            for (int k = 0; k <= pavg->nz; k++) {
                if (pavg->n_cached[i][j][k] > 0) {
                    free(pavg->particle_indices[i][j][k]);
                    free(pavg->form_factors[i][j][k]);
                }
            }
            free(pavg->n_cached[i][j]);
            free(pavg->particle_indices[i][j]);
            free(pavg->form_factors[i][j]);
        }
        free(pavg->n_cached[i]);
        free(pavg->particle_indices[i]);
        free(pavg->form_factors[i]);
    }
    free(pavg->n_cached);
    free(pavg->particle_indices);
    free(pavg->form_factors);

    /* Free main data */
    for (int i = 0; i <= pavg->nx; i++) {
        for (int j = 0; j <= pavg->ny; j++) {
            for (int k = 0; k <= pavg->nz; k++) {
                aligned_free(pavg->data[i][j][k]);
            }
            free(pavg->data[i][j]);
        }
        free(pavg->data[i]);
    }
    free(pavg->data);

    free(pavg);
}

/*
 * Clear cached particle mappings
 */
void p_avg_clear_cache(PAvg *pavg)
{
    for (int i = 0; i <= pavg->nx; i++) {
        for (int j = 0; j <= pavg->ny; j++) {
            for (int k = 0; k <= pavg->nz; k++) {
                if (pavg->n_cached[i][j][k] > 0) {
                    free(pavg->particle_indices[i][j][k]);
                    free(pavg->form_factors[i][j][k]);
                    pavg->particle_indices[i][j][k] = NULL;
                    pavg->form_factors[i][j][k] = NULL;
                    pavg->n_cached[i][j][k] = 0;
                }
            }
        }
    }
}

/*
 * Calculate local average momentum
 *
 * Direct port of p_avg_alpha() from legacy lib/p_avg_alpha.c
 */
void p_avg_calculate(PAvg *pavg, const Grid *grid, const Particle *particles,
                     int n_particles, const DrhoDr *drho_cache,
                     int center_ix, int center_iy, int center_iz)
{
    /* Mode: clear cache and return */
    if (center_ix < 0) {
        p_avg_clear_cache(pavg);
        return;
    }

    /* Determine update range */
    int i_min, i_max, j_min, j_max, k_min, k_max;

    if (center_ix == 0) {
        /* Full grid update */
        i_min = FF_RANGE + 1;
        j_min = FF_RANGE + 1;
        k_min = FF_RANGE + 1;
        i_max = pavg->nx - FF_RANGE;
        j_max = pavg->ny - FF_RANGE;
        k_max = pavg->nz - FF_RANGE;
    } else {
        /* Localized update: ±4*ff_range around center */
        int temp;

        temp = center_ix - 4 * FF_RANGE;
        i_min = (temp > FF_RANGE) ? temp : FF_RANGE + 1;

        temp = center_iy - 4 * FF_RANGE;
        j_min = (temp > FF_RANGE) ? temp : FF_RANGE + 1;

        temp = center_iz - 4 * FF_RANGE;
        k_min = (temp > FF_RANGE) ? temp : FF_RANGE + 1;

        temp = center_ix + 4 * FF_RANGE;
        i_max = (temp <= pavg->nx - FF_RANGE) ? temp : pavg->nx - FF_RANGE;

        temp = center_iy + 4 * FF_RANGE;
        j_max = (temp <= pavg->ny - FF_RANGE) ? temp : pavg->ny - FF_RANGE;

        temp = center_iz + 4 * FF_RANGE;
        k_max = (temp <= pavg->nz - FF_RANGE) ? temp : pavg->nz - FF_RANGE;
    }

    /* Build grid_tag from particles (needed for particle lookup) */
    /* This maps each grid cell to particles within it */
    int ***grid_tag_count = malloc((pavg->nx + 1) * sizeof(int**));
    int ****grid_tag = malloc((pavg->nx + 1) * sizeof(int***));

    for (int i = 0; i <= pavg->nx; i++) {
        grid_tag_count[i] = calloc(pavg->ny + 1, sizeof(int*));
        grid_tag[i] = malloc((pavg->ny + 1) * sizeof(int**));
        for (int j = 0; j <= pavg->ny; j++) {
            grid_tag_count[i][j] = calloc(pavg->nz + 1, sizeof(int));
            grid_tag[i][j] = malloc((pavg->nz + 1) * sizeof(int*));
            for (int k = 0; k <= pavg->nz; k++) {
                grid_tag[i][j][k] = malloc(PARTICLES_PER_GRID * sizeof(int));
            }
        }
    }

    /* Assign particles to grid cells */
    for (int p = 0; p < n_particles; p++) {
        int ix, iy, iz;
        grid_coords(grid, particles[p].x, particles[p].y, particles[p].z,
                    &ix, &iy, &iz);

        /* Bounds check */
        if (ix < 0 || ix > pavg->nx || iy < 0 || iy > pavg->ny ||
            iz < 0 || iz > pavg->nz) {
            continue;
        }

        /* Store particle index (using 1-based indexing for compatibility) */
        int count = grid_tag_count[ix][iy][iz];
        if (count < PARTICLES_PER_GRID - 1) {
            grid_tag[ix][iy][iz][count] = p + 1;  /* 1-based */
            grid_tag_count[ix][iy][iz]++;
        }
    }

    /* Loop over grid sites - PARALLELIZED */
    #pragma omp parallel for collapse(3) schedule(dynamic)
    for (int i = i_min; i <= i_max; i++) {
        for (int j = j_min; j <= j_max; j++) {
            for (int k = k_min; k <= k_max; k++) {

                /* Check if we already have cached data for this site */
                if (pavg->n_cached[i][j][k] == 0) {
                    /* No cache - need to compute */

                    /* Check if this site has nonzero density */
                    if (grid->baryon_density[grid_index(grid, i, j, k)] == 0.0) {
                        continue;
                    }

                    /* Temporary storage for particles contributing to this site */
                    int *temp_particles = malloc(MAX_PARTICLES * sizeof(int));
                    double *temp_g = malloc(MAX_PARTICLES * sizeof(double));
                    int n_found = 0;

                    /* Search over volume encompassing form factor range */
                    int ii_min = i - FF_RANGE;
                    int jj_min = j - FF_RANGE;
                    int kk_min = k - FF_RANGE;
                    int ii_max = i + FF_RANGE;
                    int jj_max = j + FF_RANGE;
                    int kk_max = k + FF_RANGE;

                    double px_sum = 0.0, py_sum = 0.0, pz_sum = 0.0, g_total = 0.0;

                    for (int ii = ii_min; ii <= ii_max; ii++) {
                        for (int jj = jj_min; jj <= jj_max; jj++) {
                            for (int kk = kk_min; kk <= kk_max; kk++) {
                                /* Get particles in this cell */
                                int m = grid_tag_count[ii][jj][kk];
                                for (int n = 0; n < m; n++) {
                                    int L = grid_tag[ii][jj][kk][n];  /* 1-based particle index */

                                    /* Get form factor contribution */
                                    double g = get_form_factor(drho_cache, L, i, j, k, ii, jj, kk);

                                    if (g > 0.0) {
                                        /* Store for caching */
                                        temp_particles[n_found] = L;
                                        temp_g[n_found] = g;
                                        n_found++;

                                        /* Accumulate weighted momentum */
                                        int pidx = L - 1;  /* Convert to 0-based */
                                        g_total += g;
                                        px_sum += g * particles[pidx].px;
                                        py_sum += g * particles[pidx].py;
                                        pz_sum += g * particles[pidx].pz;
                                    }
                                }
                            }
                        }
                    }

                    /* Normalize and store result */
                    if (g_total > 0.0) {
                        double norm = 1.0 / g_total;
                        pavg->data[i][j][k][0] = px_sum * norm;
                        pavg->data[i][j][k][1] = py_sum * norm;
                        pavg->data[i][j][k][2] = pz_sum * norm;
                    }

                    /* Cache the particle list for future use */
                    if (n_found > 0) {
                        pavg->n_cached[i][j][k] = n_found;
                        pavg->particle_indices[i][j][k] = malloc(n_found * sizeof(int));
                        pavg->form_factors[i][j][k] = malloc(n_found * sizeof(double));
                        memcpy(pavg->particle_indices[i][j][k], temp_particles, n_found * sizeof(int));
                        memcpy(pavg->form_factors[i][j][k], temp_g, n_found * sizeof(double));
                    }

                    free(temp_particles);
                    free(temp_g);

                } else {
                    /* Have cached data - reuse particle list */
                    int n_cached = pavg->n_cached[i][j][k];
                    double px_sum = 0.0, py_sum = 0.0, pz_sum = 0.0, g_total = 0.0;

                    for (int n = 0; n < n_cached; n++) {
                        int L = pavg->particle_indices[i][j][k][n];  /* 1-based */
                        double g = pavg->form_factors[i][j][k][n];

                        int pidx = L - 1;  /* Convert to 0-based */
                        g_total += g;
                        px_sum += g * particles[pidx].px;
                        py_sum += g * particles[pidx].py;
                        pz_sum += g * particles[pidx].pz;
                    }

                    /* Normalize and store result */
                    if (g_total > 0.0) {
                        double norm = 1.0 / g_total;
                        pavg->data[i][j][k][0] = px_sum * norm;
                        pavg->data[i][j][k][1] = py_sum * norm;
                        pavg->data[i][j][k][2] = pz_sum * norm;
                    }
                }
            }
        }
    }

    /* Free grid_tag */
    for (int i = 0; i <= pavg->nx; i++) {
        for (int j = 0; j <= pavg->ny; j++) {
            for (int k = 0; k <= pavg->nz; k++) {
                free(grid_tag[i][j][k]);
            }
            free(grid_tag_count[i][j]);
            free(grid_tag[i][j]);
        }
        free(grid_tag_count[i]);
        free(grid_tag[i]);
    }
    free(grid_tag_count);
    free(grid_tag);
}

/*
 * Get average momentum at particle position
 */
void p_avg_at_particle(const PAvg *pavg, const Grid *grid,
                       const Particle *particle,
                       double *px_avg, double *py_avg, double *pz_avg)
{
    int ix, iy, iz;
    grid_coords(grid, particle->x, particle->y, particle->z, &ix, &iy, &iz);

    /* Bounds check */
    if (ix < 0 || ix > pavg->nx || iy < 0 || iy > pavg->ny ||
        iz < 0 || iz > pavg->nz) {
        *px_avg = 0.0;
        *py_avg = 0.0;
        *pz_avg = 0.0;
        return;
    }

    /* Return local average momentum at particle's grid cell */
    *px_avg = pavg->data[ix][iy][iz][0];
    *py_avg = pavg->data[ix][iy][iz][1];
    *pz_avg = pavg->data[ix][iy][iz][2];
}
