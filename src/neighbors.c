#include "neighbors.h"
#include "grid.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

NeighborSystem* neighbors_create(size_t num_particles, float cutoff, int rebuild_interval) {
    NeighborSystem *ns = malloc(sizeof(NeighborSystem));
    if (!ns) return NULL;

    ns->num_particles = num_particles;
    ns->cutoff = cutoff;
    ns->rebuild_interval = rebuild_interval;
    ns->steps_since_rebuild = 0;

    ns->lists = calloc(num_particles, sizeof(NeighborList));
    if (!ns->lists) {
        free(ns);
        return NULL;
    }

    /* Pre-allocate space for typical neighbor count */
    int typical_neighbors = 150;  /* Conservative estimate */
    for (size_t i = 0; i < num_particles; i++) {
        ns->lists[i].capacity = typical_neighbors;
        ns->lists[i].indices = malloc(typical_neighbors * sizeof(int));
        ns->lists[i].weights = malloc(typical_neighbors * sizeof(float));
        ns->lists[i].count = 0;

        if (!ns->lists[i].indices || !ns->lists[i].weights) {
            neighbors_free(ns);
            return NULL;
        }
    }

    return ns;
}

void neighbors_free(NeighborSystem *ns) {
    if (!ns) return;

    if (ns->lists) {
        for (size_t i = 0; i < ns->num_particles; i++) {
            free(ns->lists[i].indices);
            free(ns->lists[i].weights);
        }
        free(ns->lists);
    }
    free(ns);
}

/* Form factor (Monaghan cubic spline) */
static float form_factor_1d(float r, float h) {
    float q = fabsf(r) / h;
    if (q >= 2.0f) return 0.0f;
    if (q >= 1.0f) {
        float t = 2.0f - q;
        return 0.25f * t * t * t;
    }
    return 1.0f - 1.5f * q * q + 0.75f * q * q * q;
}

static float form_factor_3d(float dx, float dy, float dz, float h) {
    return form_factor_1d(dx, h) *
           form_factor_1d(dy, h) *
           form_factor_1d(dz, h);
}

void neighbors_build(NeighborSystem *ns, const SimState *state) {
    const Grid *grid = &state->config->grid;
    const Particle *particles = state->particles;
    float h = grid->dx;  /* Form factor width */

    /* Clear existing lists */
    for (size_t i = 0; i < ns->num_particles; i++) {
        ns->lists[i].count = 0;
    }

    /* Build neighbor lists using grid cells for efficiency */
    for (size_t i = 0; i < state->num_particles; i++) {
        const Particle *pi = &particles[i];
        NeighborList *list = &ns->lists[i];

        /* Get cell coordinates */
        int ix, iy, iz;
        grid_coords(grid, pi->x, pi->y, pi->z, &ix, &iy, &iz);

        /* Search neighboring cells */
        int cell_range = (int)(ns->cutoff / grid->dx) + 1;

        for (int dx = -cell_range; dx <= cell_range; dx++) {
            for (int dy = -cell_range; dy <= cell_range; dy++) {
                for (int dz = -cell_range; dz <= cell_range; dz++) {
                    int cx = ix + dx;
                    int cy = iy + dy;
                    int cz = iz + dz;

                    /* Check bounds */
                    if (cx < 0 || cx >= grid->nx ||
                        cy < 0 || cy >= grid->ny ||
                        cz < 0 || cz >= grid->nz) {
                        continue;
                    }

                    int cell_idx = grid_index(grid, cx, cy, cz);
                    int cell_count = grid->cell_counts[cell_idx];

                    /* Check particles in this cell */
                    for (int k = 0; k < cell_count; k++) {
                        size_t j = grid->cell_particles[cell_idx][k];
                        if (j == i) continue;  /* Skip self */

                        const Particle *pj = &particles[j];

                        /* Calculate distance */
                        float rx = pj->x - pi->x;
                        float ry = pj->y - pi->y;
                        float rz = pj->z - pi->z;
                        float r = sqrtf(rx*rx + ry*ry + rz*rz);

                        if (r < ns->cutoff) {
                            /* Add to neighbor list */
                            if (list->count >= list->capacity) {
                                /* Resize */
                                list->capacity *= 2;
                                list->indices = realloc(list->indices,
                                                       list->capacity * sizeof(int));
                                list->weights = realloc(list->weights,
                                                       list->capacity * sizeof(float));
                            }

                            list->indices[list->count] = j;
                            list->weights[list->count] = form_factor_3d(rx, ry, rz, h);
                            list->count++;
                        }
                    }
                }
            }
        }
    }

    ns->steps_since_rebuild = 0;
}

int neighbors_needs_rebuild(const NeighborSystem *ns) {
    return ns->steps_since_rebuild >= ns->rebuild_interval;
}

void neighbors_step(NeighborSystem *ns) {
    ns->steps_since_rebuild++;
}
