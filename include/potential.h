#ifndef LHBUU_POTENTIAL_H
#define LHBUU_POTENTIAL_H

#include "types.h"

/* Calculate mean-field potential on grid */
void potential_calculate_meanfield(Grid *grid, const PotentialParams *params);

/* Get single-particle potential energy at grid point */
float potential_at_point(const Grid *grid, const PotentialParams *params, int ix, int iy, int iz);

/* Calculate force on particle from potentials */
void potential_calculate_force(const Grid *grid, const PotentialParams *params,
                                const Particle *p, float *fx, float *fy, float *fz);

/* Momentum-dependent potential */
float potential_momentum_dependent(const Grid *grid, const PotentialParams *params,
                                   const Particle *p, int ix, int iy, int iz);

/* Gradient of momentum-dependent potential */
void potential_grad_momentum_dep(const Grid *grid, const PotentialParams *params,
                                 const Particle *p, int ix, int iy, int iz,
                                 float *gx, float *gy, float *gz);

#endif /* LHBUU_POTENTIAL_H */
