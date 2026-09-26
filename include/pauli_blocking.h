#ifndef LHBUU_PAULI_BLOCKING_H
#define LHBUU_PAULI_BLOCKING_H

#include "types.h"

/* Pauli blocking constants */
#define DEGENERACY_SPIN_ISOSPIN 4  /* g=4 for spin and isospin */
#define DEGENERACY_SPIN 2          /* g=2 for spin only */

/* Pauli blocking state (initialized once per simulation) */
typedef struct {
    float r_pauli;              /* Configuration space Pauli radius (fm) */
    float p_pauli;              /* Momentum space Pauli radius (GeV) */
    float r_pauli2;             /* r_pauli squared */
    float p_pauli2;             /* p_pauli squared */
    float max_occupancy;        /* Maximum phase space occupancy */

    /* Isospin-dependent factors */
    float r_neutron_iso_factor; /* Neutron config space multiplier */
    float p_neutron_iso_factor; /* Neutron momentum space multiplier */
    float r_proton_iso_factor;  /* Proton config space multiplier */
    float p_proton_iso_factor;  /* Proton momentum space multiplier */

    /* Cold nucleus tracking */
    int *proj_targ_interaction; /* Array: has particle hit opposite nucleus? */
} PauliState;

/**
 * Initialize Pauli blocking radii and parameters.
 *
 * Physics: Pauli blocking prevents collisions into occupied phase space cells.
 * Uses 6D phase space (3 position + 3 momentum) with radii determined by
 * Fermi momentum and nuclear radius.
 *
 * @param state Simulation state
 * @param n_pauli_ensembles Number of ensembles to sample for blocking check
 * @param n_full Full ensemble count for phase space normalization
 * @param pf0 Fermi momentum at saturation density (GeV)
 * @param diff_rn_rp Neutron-proton radius difference (fm)
 * @return Initialized PauliState or NULL on error
 */
PauliState* pauli_blocking_init(const SimState *state,
                                 int n_pauli_ensembles,
                                 int n_full,
                                 float pf0,
                                 float diff_rn_rp);

/**
 * Free Pauli blocking state.
 */
void pauli_blocking_free(PauliState *ps);

/**
 * Check if collision is Pauli blocked.
 *
 * Algorithm:
 * 1. Apply cold nucleus filter: reject collisions within same cold nucleus
 * 2. For each particle in collision pair:
 *    - Sample random ensembles (or all if n_pauli_ensembles == N_ensembles)
 *    - Use grid to find nearby particles within r_pauli
 *    - Count particles in same phase space cell (within r_pauli and p_pauli)
 *    - Apply isospin-dependent radii if D_pot != 0
 * 3. Accept/reject based on occupancy vs random number
 *
 * @param state Simulation state
 * @param ps Pauli blocking state
 * @param p1 First particle
 * @param p2 Second particle
 * @return 1 if blocked, 0 if allowed
 */
int pauli_blocking_check(const SimState *state,
                          PauliState *ps,
                          const Particle *p1,
                          const Particle *p2);

/**
 * Mark particle as having interacted with opposite nucleus.
 * Used for cold nucleus filter.
 */
void pauli_blocking_mark_interaction(PauliState *ps, int particle_idx);

#endif /* LHBUU_PAULI_BLOCKING_H */
