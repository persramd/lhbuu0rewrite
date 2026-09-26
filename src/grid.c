#include "grid.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int grid_init(Grid *grid, int nx, int ny, int nz, float dx) {
    grid->nx = nx;
    grid->ny = ny;
    grid->nz = nz;
    grid->dx = dx;

    /* Center grid at origin */
    grid->origin_x = -0.5f * nx * dx;
    grid->origin_y = -0.5f * ny * dx;
    grid->origin_z = -0.5f * nz * dx;

    size_t n_cells = nx * ny * nz;

    grid->baryon_density = calloc(n_cells, sizeof(float));
    grid->charge_density = calloc(n_cells, sizeof(float));
    grid->potential = calloc(n_cells, sizeof(float));
    grid->coulomb_potential = calloc(n_cells, sizeof(float));

    if (!grid->baryon_density || !grid->charge_density ||
        !grid->potential || !grid->coulomb_potential) {
        error_exit("Failed to allocate grid arrays");
    }

    grid->max_per_cell = 100;  /* Initial estimate */
    grid->cell_particles = malloc(n_cells * sizeof(int*));
    grid->cell_counts = calloc(n_cells, sizeof(int));

    for (size_t i = 0; i < n_cells; i++) {
        grid->cell_particles[i] = malloc(grid->max_per_cell * sizeof(int));
    }

    return 1;
}

void grid_free(Grid *grid) {
    free(grid->baryon_density);
    free(grid->charge_density);
    free(grid->potential);
    free(grid->coulomb_potential);
    free(grid->cell_counts);

    size_t n_cells = grid->nx * grid->ny * grid->nz;
    for (size_t i = 0; i < n_cells; i++) {
        free(grid->cell_particles[i]);
    }
    free(grid->cell_particles);
}

/* Form factor (Monaghan cubic spline kernel) */
static float form_factor_1d(float r) {
    float q = fabsf(r);
    if (q >= 2.0f) return 0.0f;
    if (q >= 1.0f) return 0.25f * (2.0f - q) * (2.0f - q) * (2.0f - q);
    return 1.0f - 1.5f * q * q + 0.75f * q * q * q;
}

void grid_calculate_density(Grid *grid, const Particle *particles,
                             size_t n_particles, int ff_range) {
    /* Zero density grids */
    size_t n_cells = grid->nx * grid->ny * grid->nz;
    memset(grid->baryon_density, 0, n_cells * sizeof(float));
    memset(grid->charge_density, 0, n_cells * sizeof(float));

    /* Accumulate density from particles (parallelized) */
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < n_particles; i++) {
        const Particle *p = &particles[i];

        int ix, iy, iz;
        grid_coords(grid, p->x, p->y, p->z, &ix, &iy, &iz);

        /* Distribute over neighboring cells */
        for (int dx = -ff_range; dx <= ff_range; dx++) {
            for (int dy = -ff_range; dy <= ff_range; dy++) {
                for (int dz = -ff_range; dz <= ff_range; dz++) {
                    int cx = ix + dx;
                    int cy = iy + dy;
                    int cz = iz + dz;

                    /* Periodic boundary conditions */
                    if (cx < 0 || cx >= grid->nx) continue;
                    if (cy < 0 || cy >= grid->ny) continue;
                    if (cz < 0 || cz >= grid->nz) continue;

                    float ff_x = form_factor_1d((float)dx);
                    float ff_y = form_factor_1d((float)dy);
                    float ff_z = form_factor_1d((float)dz);
                    float ff = ff_x * ff_y * ff_z;

                    int idx = grid_index(grid, cx, cy, cz);

                    /* Atomic update for thread safety */
                    #pragma omp atomic
                    grid->baryon_density[idx] += ff;
                    if (p->charge) {
                        #pragma omp atomic
                        grid->charge_density[idx] += ff;
                    }
                }
            }
        }
    }

    /* Normalize by cell volume */
    float cell_volume = grid->dx * grid->dx * grid->dx;
    for (size_t i = 0; i < n_cells; i++) {
        grid->baryon_density[i] /= cell_volume;
        grid->charge_density[i] /= cell_volume;
    }
}

void grid_assign_particles(Grid *grid, Particle *particles, size_t n_particles) {
    size_t n_cells = grid->nx * grid->ny * grid->nz;

    /* Clear counts */
    memset(grid->cell_counts, 0, n_cells * sizeof(int));

    /* Assign particles to cells */
    for (size_t i = 0; i < n_particles; i++) {
        Particle *p = &particles[i];

        int ix, iy, iz;
        grid_coords(grid, p->x, p->y, p->z, &ix, &iy, &iz);

        if (ix < 0 || ix >= grid->nx ||
            iy < 0 || iy >= grid->ny ||
            iz < 0 || iz >= grid->nz) {
            continue;  /* Particle outside grid */
        }

        int idx = grid_index(grid, ix, iy, iz);
        int count = grid->cell_counts[idx];

        if (count >= grid->max_per_cell) {
            warning("Cell overflow at (%d,%d,%d)", ix, iy, iz);
            continue;
        }

        grid->cell_particles[idx][count] = i;
        grid->cell_counts[idx]++;
        p->grid_cell = idx;
    }
}

void grid_solve_coulomb(Grid *grid) {
    /* Simplified Coulomb solver - direct sum (slow but simple) */
    /* For production, use FFT-based Poisson solver */

    size_t n_cells = grid->nx * grid->ny * grid->nz;
    memset(grid->coulomb_potential, 0, n_cells * sizeof(float));

    float alpha = 1.44f;  /* e^2 in MeV·fm */

    for (int ix = 0; ix < grid->nx; ix++) {
        for (int iy = 0; iy < grid->ny; iy++) {
            for (int iz = 0; iz < grid->nz; iz++) {
                int idx = grid_index(grid, ix, iy, iz);
                float phi = 0.0f;

                /* Sum over all cells */
                for (int jx = 0; jx < grid->nx; jx++) {
                    for (int jy = 0; jy < grid->ny; jy++) {
                        for (int jz = 0; jz < grid->nz; jz++) {
                            if (ix == jx && iy == jy && iz == jz) continue;

                            int jdx = grid_index(grid, jx, jy, jz);
                            float rho = grid->charge_density[jdx];
                            if (rho == 0.0f) continue;

                            float dx = (ix - jx) * grid->dx;
                            float dy = (iy - jy) * grid->dx;
                            float dz = (iz - jz) * grid->dx;
                            float r = sqrtf(dx*dx + dy*dy + dz*dz);

                            if (r > 0.0f) {
                                float dV = grid->dx * grid->dx * grid->dx;
                                phi += alpha * rho * dV / r;
                            }
                        }
                    }
                }

                grid->coulomb_potential[idx] = phi;
            }
        }
    }
}
