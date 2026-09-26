#ifndef LHBUU_INTEGRATOR_H
#define LHBUU_INTEGRATOR_H

#include "types.h"

/* Velocity Verlet integration step */
void integrator_step(SimState *state);

/* Update positions (first half of Verlet) */
void integrator_update_positions(SimState *state);

/* Update momenta (second half of Verlet) */
void integrator_update_momenta(SimState *state);

/* Calculate forces on all particles */
void integrator_calculate_forces(SimState *state);

#endif /* LHBUU_INTEGRATOR_H */
