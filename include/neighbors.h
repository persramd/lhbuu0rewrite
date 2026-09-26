#ifndef LHBUU_NEIGHBORS_H
#define LHBUU_NEIGHBORS_H

#include "types.h"

/* Neighbor list for efficient local calculations */
typedef struct {
    int *indices;           /* Neighbor particle indices */
    float *weights;         /* Pre-computed form factor weights */
    int count;              /* Number of neighbors */
    int capacity;           /* Allocated capacity */
} NeighborList;

/* Neighbor list collection for all particles */
typedef struct {
    NeighborList *lists;    /* One list per particle */
    size_t num_particles;
    float cutoff;           /* Cutoff radius */
    int rebuild_interval;   /* Timesteps between rebuilds */
    int steps_since_rebuild;
} NeighborSystem;

/* Create neighbor list system */
NeighborSystem* neighbors_create(size_t num_particles, float cutoff, int rebuild_interval);

/* Free neighbor list system */
void neighbors_free(NeighborSystem *ns);

/* Build/rebuild neighbor lists */
void neighbors_build(NeighborSystem *ns, const SimState *state);

/* Check if rebuild is needed */
int neighbors_needs_rebuild(const NeighborSystem *ns);

/* Update rebuild counter */
void neighbors_step(NeighborSystem *ns);

#endif /* LHBUU_NEIGHBORS_H */
