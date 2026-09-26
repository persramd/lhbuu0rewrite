#ifndef LHBUU_INTEGRATOR_H
#define LHBUU_INTEGRATOR_H

#include "types.h"

/* Modified Velocity Verlet integration step
 *
 * Algorithm (Allen & Tildesley):
 *   1. Calculate F(t) = -∇_r V at current positions
 *   2. Half-step momentum: p_half = p(t) + 0.5*F(t)*dt
 *   3. Calculate ∇_p V(r,p) for momentum-dependent potential
 *   4. Full-step position: r(t+dt) = r(t) + (p_half/E)*dt + ∇_p V*dt
 *   5. Update densities and potentials at r(t+dt)
 *   6. Calculate F(t+dt) = -∇_r V at new positions
 *   7. Half-step momentum corrector: p(t+dt) = p_half + 0.5*F(t+dt)*dt
 *
 * Key: Position update includes momentum-dependent gradient term
 */
void integrator_step(SimState *state);

/* Calculate forces on all particles: F = -∇_r V */
void integrator_calculate_forces(SimState *state);

/* Free integrator temporary storage */
void integrator_cleanup(void);

#endif /* LHBUU_INTEGRATOR_H */
