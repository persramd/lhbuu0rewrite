/*
 * inside_sigma.c - Collision geometric pre-filter
 *
 * Checks all pairs of nucleons from a given ensemble and decides
 * whether a collision is possible. The test is passed if:
 *
 * 1) Two nucleons are separated by a distance less than the distance
 *    corresponding to the nucleon-nucleon cross section (takes into
 *    account the maximum possible distance a nucleon can travel in
 *    one time step)
 *
 * 2) The nucleons pass the distance of closest approach in the current
 *    time step. To calculate this, we assume that the momenta of the
 *    nucleons in question is the same at the start and end of the
 *    time step.
 *
 * FUNCTION ONLY FOR nucleon-nucleon ELASTIC SCATTERING
 *
 * Note: This function does not check that the (Lorentz contracted)
 * closest approach in the nucleon-nucleon collision CM frame is within
 * the cross section. (see collision.c for this)
 *
 * This function is a "first pass test" for scattering test. It uses
 * the MAXIMUM total cross section to test for an allowed collision.
 * A finer test is performed in function/file 'scatter.c'
 *
 * Physics preserved from original LHBUU (1997)
 */

#include <math.h>
#include "inside_sigma.h"
#include "constants.h"
#include "macros.h"

/*
 * Check if two particles can possibly collide geometrically
 *
 * Algorithm (lines 80-103 from original):
 * 1. Box test: Check if |Δx|, |Δy|, |Δz| < delta_r_max
 *    - Uses box instead of sphere for speed (overcounts slightly)
 *    - Correct sphere check done later in collision.c
 *
 * 2. Closest approach test:
 *    - Compute r·p at start and end of timestep
 *    - If signs differ, particles passed closest approach
 *    - If r·p_start * r·p_end < 0, collision possible
 *
 * delta_r_max = sqrt(sigma_max/(π*10) + 4*dt²)
 *   - First term: radius from cross section (sigma in mb, factor 10 for units)
 *   - Second term: 4*dt² accounts for maximum relative motion
 */
int inside_sigma(const Particle *p1, const Particle *p2,
                 double x_new1, double y_new1, double z_new1,
                 double x_new2, double y_new2, double z_new2,
                 double sigma_nn_max, double dt)
{
    /* Separation distance for variable cross section
     * Original: delta_r_max_nucnuc = sqrtf(sigma_nn_max_local/(PIE*10) + 4*square(dt))
     * Factor of 10: converts mb to fm² (1 mb = 0.1 fm²)
     */
    double delta_r_max_nucnuc = sqrt(sigma_nn_max / (M_PI * 10.0) + 4.0 * square(dt));

    /* Configuration space separations at start of timestep
     * NOTE: This method overcounts as particles are checked for inside the
     * cross section by using a box, so the edges of the box are too big.
     * However, the correct sphere checking is done in collisions.c
     */

    /* Box test in x-direction */
    double delta_rx_start = p1->x - p2->x;
    if (fabs(delta_rx_start) > delta_r_max_nucnuc) {
        return 0;
    }

    /* Box test in y-direction */
    double delta_ry_start = p1->y - p2->y;
    if (fabs(delta_ry_start) > delta_r_max_nucnuc) {
        return 0;
    }

    /* Box test in z-direction */
    double delta_rz_start = p1->z - p2->z;
    if (fabs(delta_rz_start) > delta_r_max_nucnuc) {
        return 0;
    }

    /* Check that particles pass closest approach. If they do, the sign of the
     * dot product (delta_r, delta_p) changes sign
     *
     * r·p at start of timestep
     */
    double delta_px = p1->px - p2->px;
    double delta_py = p1->py - p2->py;
    double delta_pz = p1->pz - p2->pz;

    double r_dot_p_start = delta_rx_start * delta_px
                         + delta_ry_start * delta_py
                         + delta_rz_start * delta_pz;

    /* r·p at end of timestep */
    double r_dot_p_end = (x_new1 - x_new2) * delta_px
                       + (y_new1 - y_new2) * delta_py
                       + (z_new1 - z_new2) * delta_pz;

    /* If r·p changes sign, particles passed closest approach
     * Original: if(r_dot_p_start*r_dot_p_end<0) { possible_collisions[i][j]=1; }
     */
    if (r_dot_p_start * r_dot_p_end < 0.0) {
        return 1;  /* Collision possible */
    } else {
        return 0;  /* Collision impossible */
    }
}
