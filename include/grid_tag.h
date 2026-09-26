/*
 * grid_tag: Spatial Hashing for Fast Neighbor Lookup
 *
 * Maps particles to grid cells for O(1) neighbor search
 * Critical for momentum-dependent potentials and collision detection
 *
 * From original lib/grid_sort.c
 */

#ifndef GRID_TAG_H
#define GRID_TAG_H

#include "constants.h"

/*
 * GridTag structure
 *
 * 4D array: grid_tag[i][j][k][particle_list]
 * - grid_tag[i][j][k][0] = number of particles in cell (i,j,k)
 * - grid_tag[i][j][k][1..count] = particle indices in this cell
 *
 * Memory: NX × NY × NZ × PARTICLES_PER_GRID × 4 bytes
 *         64 × 64 × 64 × 16 × 4 = 16 MB (acceptable)
 */
typedef struct {
    int ****data;           /* [i][j][k][particle_list] */
    int dims[3];            /* Grid dimensions [NX, NY, NZ] */
    int max_particles;      /* Max particles per cell */
    int max_count;          /* Actual max particles seen in any cell */
} GridTag;

/*
 * Initialize grid_tag structure
 *
 * Args:
 *   nx, ny, nz: grid dimensions
 *   max_particles_per_cell: allocation size per cell
 *
 * Returns:
 *   Allocated GridTag structure
 */
GridTag* grid_tag_init(int nx, int ny, int nz, int max_particles_per_cell);

/*
 * Free grid_tag structure
 */
void grid_tag_free(GridTag *grid);

/*
 * Assign particles to grid cells
 *
 * Rebuilds cell membership based on current particle positions
 * Called every timestep after position update
 *
 * Args:
 *   grid: grid_tag structure
 *   n_particles: total number of particles
 *   x, y, z: particle positions [1..n_particles] (1-based indexing)
 *   grid_center_x, grid_center_y, grid_center_z: grid origins
 *
 * From original: lib/grid_sort.c lines 21-68 (assign_grid_tags)
 *
 * Note: Does NOT physically reorder particles (sort_count_def=0 in production)
 */
void grid_tag_assign(GridTag *grid, int n_particles,
                     const double *x, const double *y, const double *z,
                     double grid_center_x, double grid_center_y, double grid_center_z);

/*
 * Clear grid_tag (zero all counts)
 *
 * Called before rebuilding
 */
void grid_tag_clear(GridTag *grid);

/*
 * Get number of particles in cell
 *
 * Args:
 *   grid: grid_tag structure
 *   i, j, k: grid cell indices
 *
 * Returns:
 *   Number of particles in cell (i,j,k)
 */
static inline int grid_tag_count(const GridTag *grid, int i, int j, int k)
{
    return grid->data[i][j][k][0];
}

/*
 * Get particle index in cell
 *
 * Args:
 *   grid: grid_tag structure
 *   i, j, k: grid cell indices
 *   n: particle number in list (1..count)
 *
 * Returns:
 *   Particle index
 */
static inline int grid_tag_particle(const GridTag *grid, int i, int j, int k, int n)
{
    return grid->data[i][j][k][n];
}

#endif /* GRID_TAG_H */
