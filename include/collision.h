#ifndef LHBUU_COLLISION_H
#define LHBUU_COLLISION_H

#include "types.h"

/* Process all collisions for current timestep */
int collision_process_all(SimState *state);

/* Check if two particles can collide (geometric test) */
int collision_check_geometric(const Particle *p1, const Particle *p2, float sigma_max);

/* Perform scattering calculation */
int collision_scatter(Particle *p1, Particle *p2, const CollisionParams *params);

/* Check Pauli blocking for final state */
int collision_check_pauli(const SimState *state, const Particle *p1, const Particle *p2);

/* Calculate NN cross section */
float collision_cross_section(const Particle *p1, const Particle *p2, const CollisionParams *params);

#endif /* LHBUU_COLLISION_H */
