#ifndef LHBUU_GRID_H
#define LHBUU_GRID_H

#include "types.h"

/* Grid initialization and cleanup */
int grid_init(Grid *grid, int nx, int ny, int nz, float dx);
void grid_free(Grid *grid);

/* Density calculations */
void grid_calculate_density(Grid *grid, const Particle *particles, size_t n_particles, int form_factor_range);

/* Grid indexing */
static inline int grid_index(const Grid *grid, int ix, int iy, int iz) {
    return ix * grid->ny * grid->nz + iy * grid->nz + iz;
}

static inline void grid_coords(const Grid *grid, float x, float y, float z, int *ix, int *iy, int *iz) {
    *ix = (int)((x - grid->origin_x) / grid->dx + 0.5);
    *iy = (int)((y - grid->origin_y) / grid->dx + 0.5);
    *iz = (int)((z - grid->origin_z) / grid->dx + 0.5);
}

/* Particle-cell associations */
void grid_assign_particles(Grid *grid, Particle *particles, size_t n_particles);
void grid_sort_particles(Grid *grid, Particle *particles, size_t n_particles);

/* Coulomb solver */
void grid_solve_coulomb(Grid *grid);

#endif /* LHBUU_GRID_H */
