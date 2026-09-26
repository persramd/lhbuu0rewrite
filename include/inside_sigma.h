/*
 * inside_sigma.h - Collision geometric pre-filter
 *
 * Checks if two nucleons from a given ensemble can possibly collide.
 * This is a "first pass test" using the MAXIMUM total cross section.
 *
 * Test criteria:
 * 1) Particles separated by less than delta_r_max (includes max travel distance)
 * 2) Particles pass distance of closest approach in current time step
 *
 * Physics preserved from original LHBUU (1997):
 * - Box test for computational efficiency
 * - r·p sign change test for closest approach
 * - Assumes constant momentum during time step
 *
 * Note: Does NOT check Lorentz-contracted distance in CM frame
 *       (that's done in collision_scatter.c for finer testing)
 */

#ifndef INSIDE_SIGMA_H
#define INSIDE_SIGMA_H

#include "types.h"

/*
 * Check if two particles can possibly collide geometrically
 *
 * Arguments:
 *   p1, p2: Current particle states
 *   x_new, y_new, z_new: Proposed positions at end of timestep
 *   sigma_nn_max: Maximum NN cross section (mb)
 *   dt: Time step (fm/c)
 *
 * Returns:
 *   1 if collision is possible (passed both tests)
 *   0 if collision is impossible (failed geometric pre-filter)
 */
int inside_sigma(const Particle *p1, const Particle *p2,
                 double x_new1, double y_new1, double z_new1,
                 double x_new2, double y_new2, double z_new2,
                 double sigma_nn_max, double dt);

#endif /* INSIDE_SIGMA_H */
