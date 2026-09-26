/*
 * LHBUU Medium Corrections
 *
 * In-medium modifications to collision rates for nuclear transport.
 * These corrections account for:
 * 1. Effective mass modifications (relative velocity ratio)
 * 2. Phase space availability (density of states ratio)
 *
 * Physics:
 * - Collision rate ~ σ * v_rel * ρ_states
 * - In medium: v_rel and ρ_states modified by mean field
 * - Preserves detailed balance and energy-momentum conservation
 */

#ifndef LHBUU_MEDIUM_CORRECTIONS_H
#define LHBUU_MEDIUM_CORRECTIONS_H

#include "types.h"

/*
 * Relative velocity ratio: v_free / v_medium
 *
 * Accounts for effective mass modifications from momentum-dependent potential.
 * In free space: v_rel = 2|p_cm|/m
 * In medium:     v_rel includes gradient of mean field Hamiltonian
 *
 * Formula:
 *   v_medium = |2p_cm/m + grad_p1(H) - grad_p2(H)|
 *   ratio = v_free / v_medium
 *
 * Arguments:
 *   grid      - Grid with density and potential information
 *   params    - Potential parameters for gradient calculation
 *   p1        - First particle
 *   p2        - Second particle
 *   ix1..iz1  - Grid location of particle 1
 *   ix2..iz2  - Grid location of particle 2
 *
 * Returns:
 *   Ratio of free-space to in-medium relative velocity
 */
double medium_relative_velocity_ratio(
    const Grid *grid,
    const PotentialParams *params,
    const Particle *p1,
    const Particle *p2,
    int ix1, int iy1, int iz1,
    int ix2, int iy2, int iz2
);

/*
 * Density of states ratio: 1 / |1 + correction|
 *
 * Phase space availability modification from mean field.
 * Accounts for Pauli blocking enhancement/suppression in dense matter.
 *
 * Formula:
 *   multiplier = m / (2|p_cm|)
 *   correction = multiplier * (∂U_1/∂p_1 + ∂U_2/∂p_2)
 *   ratio = 1 / |1 + correction|
 *
 * Arguments:
 *   grid            - Grid with density and potential information
 *   params          - Potential parameters
 *   p1, p2          - Particles (current momenta)
 *   ix1..iz1        - Grid location of particle 1
 *   ix2..iz2        - Grid location of particle 2
 *   px_cm, py_cm, pz_cm - Center-of-mass momentum (post-collision)
 *
 * Returns:
 *   Phase space availability ratio
 */
double medium_density_states_ratio(
    const Grid *grid,
    const PotentialParams *params,
    const Particle *p1,
    const Particle *p2,
    int ix1, int iy1, int iz1,
    int ix2, int iy2, int iz2,
    double px_cm, double py_cm, double pz_cm
);

/*
 * Helper: Gradient of momentum-dependent potential w.r.t. particle momentum
 *
 * Calculates ∇_p U(r,p) at particle location.
 * This is the velocity correction in Hamilton's equations.
 *
 * Arguments:
 *   grid       - Grid with density information
 *   params     - Potential parameters (C, lambda)
 *   p          - Particle (position and momentum)
 *   ix, iy, iz - Grid cell indices
 *   grad_px, grad_py, grad_pz - Output gradient components
 */
void medium_grad_potential_momentum(
    const Grid *grid,
    const PotentialParams *params,
    const Particle *p,
    int ix, int iy, int iz,
    double *grad_px, double *grad_py, double *grad_pz
);

/*
 * Helper: Partial derivative of lattice potential
 *
 * Calculates ∂U_lattice/∂p for a particle at given momentum.
 * Used in density of states correction.
 *
 * Arguments:
 *   grid       - Grid with density information
 *   params     - Potential parameters
 *   p          - Particle
 *   ix, iy, iz - Grid cell indices
 *   px, py, pz - Momentum at which to evaluate derivative
 *
 * Returns:
 *   Partial derivative magnitude
 */
double medium_partial_potential_momentum(
    const Grid *grid,
    const PotentialParams *params,
    const Particle *p,
    int ix, int iy, int iz,
    double px, double py, double pz
);

#endif /* LHBUU_MEDIUM_CORRECTIONS_H */
