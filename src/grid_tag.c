/*
 * grid_tag: Spatial Hashing Implementation
 *
 * Direct port from lib/grid_sort.c (assign_grid_tags only)
 * Skips sort_grid_tags (disabled in production: sort_count_def=0)
 */

#include "grid_tag.h"
#include "constants.h"
#include "macros.h"
#include <stdlib.h>
#include <stdio.h>

/*
 * Allocate 4D array for grid_tag
 */
GridTag* grid_tag_init(int nx, int ny, int nz, int max_particles_per_cell)
{
    GridTag *grid = malloc(sizeof(GridTag));

    grid->dims[0] = nx;
    grid->dims[1] = ny;
    grid->dims[2] = nz;
    grid->max_particles = max_particles_per_cell;
    grid->max_count = 0;

    /* Allocate 4D array [i][j][k][particle_list]
     * Using 1-based indexing to match original (allocate nx+1, etc.) */
    grid->data = malloc((nx + 1) * sizeof(int***));

    for (int i = 0; i <= nx; i++) {
        grid->data[i] = malloc((ny + 1) * sizeof(int**));

        for (int j = 0; j <= ny; j++) {
            grid->data[i][j] = malloc((nz + 1) * sizeof(int*));

            for (int k = 0; k <= nz; k++) {
                /* +1 for count at index 0, max_particles for particle indices */
                grid->data[i][j][k] = calloc(max_particles_per_cell + 1, sizeof(int));
            }
        }
    }

    return grid;
}

/*
 * Free grid_tag structure
 */
void grid_tag_free(GridTag *grid)
{
    if (!grid) return;

    for (int i = 0; i <= grid->dims[0]; i++) {
        for (int j = 0; j <= grid->dims[1]; j++) {
            for (int k = 0; k <= grid->dims[2]; k++) {
                free(grid->data[i][j][k]);
            }
            free(grid->data[i][j]);
        }
        free(grid->data[i]);
    }
    free(grid->data);
    free(grid);
}

/*
 * Clear grid_tag counts
 */
void grid_tag_clear(GridTag *grid)
{
    for (int i = 1; i <= grid->dims[0]; i++) {
        for (int j = 1; j <= grid->dims[1]; j++) {
            for (int k = 1; k <= grid->dims[2]; k++) {
                grid->data[i][j][k][0] = 0;
            }
        }
    }
    grid->max_count = 0;
}

/*
 * Assign particles to grid cells
 *
 * Direct port from lib/grid_sort.c lines 21-68
 */
void grid_tag_assign(GridTag *grid, int n_particles,
                     const double *x, const double *y, const double *z,
                     double grid_center_x, double grid_center_y, double grid_center_z)
{
    const int max_x_grid = grid->dims[0];
    const int max_y_grid = grid->dims[1];
    const int max_z_grid = grid->dims[2];

    /* Clear all counts (line 43-46) */
    grid_tag_clear(grid);

    /* Loop over all particles and assign to grid cells (lines 50-67) */
    for (int n = 1; n <= n_particles; n++) {
        /* Find particle's grid cell using where_am_I macro (lines 51-53) */
        int i = where_am_I(x[n], grid_center_x, INV_DX);
        int j = where_am_I(y[n], grid_center_y, INV_DX);
        int k = where_am_I(z[n], grid_center_z, INV_DX);

        /* Check if particle is on grid (lines 54-56) */
        if (i > 0 && i <= max_x_grid &&
            j > 0 && j <= max_y_grid &&
            k > 0 && k <= max_z_grid)
        {
            /* Increment particle count in this cell (line 57) */
            grid->data[i][j][k][0]++;
            int m = grid->data[i][j][k][0];

            /* Check for overflow (lines 59-64) */
            if (m > grid->max_particles) {
                fprintf(stderr, "ERROR: grid_tag overflow at cell [%d][%d][%d]\n", i, j, k);
                fprintf(stderr, "  Particle %d causes count=%d > max=%d\n",
                        n, m, grid->max_particles);
                fprintf(stderr, "  Increase PARTICLES_PER_GRID in constants.h\n");
                exit(1);
            }

            /* Track maximum count seen (line 65) */
            if (m > grid->max_count) {
                grid->max_count = m;
            }

            /* Store particle index (line 66) */
            grid->data[i][j][k][m] = n;
        }
    }
}
