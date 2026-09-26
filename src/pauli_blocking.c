/* NOTE: COLLISIONS OFF THE GRID ARE NOT BLOCKED! */
#include "pauli_blocking.h"
#include "grid.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Legacy external constant - saturation density */
#define RHO0 0.16  /* fm^-3 */

/* Random number generator (0-1) - placeholder for actual RNG */
static float random_number(int seed) {
    /* TODO: Replace with proper RNG from random_number.c */
    return (float)rand() / (float)RAND_MAX;
}

/* Legacy grid coordinate finder */
static inline int where_am_I(float pos, float origin, float dx_inv) {
    return (int)((pos - origin) * dx_inv + 0.5) + 1;
}

/* Utility: square function */
static inline float square(float x) {
    return x * x;
}

PauliState* pauli_blocking_init(const SimState *state,
                                 int n_pauli_ensembles,
                                 int n_full,
                                 float pf0,
                                 float diff_rn_rp) {
    PauliState *ps = (PauliState*)malloc(sizeof(PauliState));
    if (!ps) return NULL;

    const SimConfig *cfg = state->config;
    const Nucleus *proj = &cfg->projectile;
    const Nucleus *targ = &cfg->target;
    float D_pot = cfg->potential.D;  /* Isospin potential strength */

    /* Allocate cold nucleus interaction tracker */
    ps->proj_targ_interaction = (int*)calloc(state->max_particles, sizeof(int));
    if (!ps->proj_targ_interaction) {
        free(ps);
        return NULL;
    }

    /* Determine smaller nucleus for Pauli parameters */
    float r_small;
    int A_proj = proj->A;
    int A_targ = targ->A;

    if (A_targ < A_proj || A_proj == 0) {
        r_small = powf(3.0 / (4.0 * PI * RHO0) * A_targ, 0.333333f);
    } else {
        r_small = powf(3.0 / (4.0 * PI * RHO0) * A_proj, 0.333333f);
    }

    /* Phase space factor */
    float phase_factor = powf((float)n_full / n_pauli_ensembles
                              / square(4.0 * PI / 3.0) / DEGENERACY_SPIN_ISOSPIN,
                              0.166667f);

    /* Configuration and momentum space Pauli radii */
    ps->r_pauli2 = phase_factor * sqrtf(2.0 * PI * HBARC * r_small / pf0);
    ps->p_pauli2 = phase_factor * sqrtf(2.0 * PI * HBARC * pf0 / r_small);

    /* Maximum occupancy per phase space cell */
    ps->max_occupancy = DEGENERACY_SPIN_ISOSPIN * square(4.0 * PI / 3.0)
                        * n_pauli_ensembles
                        * powf(ps->r_pauli2 * ps->p_pauli2 / (2.0 * PI * HBARC), 3.0)
                        - 1.0;

    /* Enforce minimum occupancy */
    if (fabsf(ps->max_occupancy - 1.0) < 1.0e-6f) {
        ps->max_occupancy = 1.0;
    }

    if (ps->max_occupancy < 1.0) {
        fprintf(stderr, "ERROR: Pauli blocking sphere too small. Increase n_full.\n");
        free(ps->proj_targ_interaction);
        free(ps);
        return NULL;
    }

    ps->r_pauli = ps->r_pauli2;
    ps->p_pauli = ps->p_pauli2;
    ps->r_pauli2 *= ps->r_pauli;
    ps->p_pauli2 *= ps->p_pauli;

    /* Isospin-dependent blocking (if D_pot != 0) */
    ps->r_neutron_iso_factor = 1.0;
    ps->p_neutron_iso_factor = 1.0;
    ps->r_proton_iso_factor = 1.0;
    ps->p_proton_iso_factor = 1.0;

    if (D_pot != 0.0) {
        /* Use smaller nucleus to set parameters */
        int avg_N, avg_Z;
        if (A_proj < A_targ && A_proj != 0) {
            avg_N = A_proj;
            avg_Z = proj->Z;
        } else {
            avg_N = A_targ;
            avg_Z = targ->Z;
        }

        /* Check for zero neutron or proton numbers */
        if (avg_N == 0) avg_N = avg_Z;
        if (avg_Z == 0) avg_Z = avg_N;

        /* Neutron and proton radii */
        float r_N = r_small + diff_rn_rp / 2.0;
        float r_P = r_small - diff_rn_rp / 2.0;
        float vol_N = (4.0 * PI / 3.0) * powf(r_N, 3.0);
        float vol_P = (4.0 * PI / 3.0) * powf(r_P, 3.0);
        float rho_N = (avg_N - avg_Z) / vol_N;
        float rho_P = avg_Z / vol_P;

        /* Check for zero densities */
        if (rho_N == 0.0) rho_N = rho_P;
        if (rho_P == 0.0) rho_P = rho_N;

        /* Fermi momenta */
        float pfermi_N = powf(6.0 * square(PI) / DEGENERACY_SPIN * rho_N, 0.333333f) * HBARC;
        float pfermi_P = powf(6.0 * square(PI) / DEGENERACY_SPIN * rho_P, 0.333333f) * HBARC;

        /* Neutron Pauli radii */
        float rpauli_N = phase_factor
                         * powf((float)DEGENERACY_SPIN_ISOSPIN / DEGENERACY_SPIN, 0.166667f)
                         * sqrtf(2.0 * PI * HBARC * r_N / pfermi_N);
        ps->r_neutron_iso_factor = rpauli_N / ps->r_pauli;

        float ppauli_N = rpauli_N * pfermi_N / r_N;
        ps->p_neutron_iso_factor = ppauli_N / ps->p_pauli;

        /* Proton Pauli radii */
        float rpauli_P = phase_factor
                         * powf((float)DEGENERACY_SPIN_ISOSPIN / DEGENERACY_SPIN, 0.166667f)
                         * sqrtf(2.0 * PI * HBARC * r_P / pfermi_P);
        ps->r_proton_iso_factor = rpauli_P / ps->r_pauli;

        float ppauli_P = rpauli_P * pfermi_P / r_P;
        ps->p_proton_iso_factor = ppauli_P / ps->p_pauli;
    }

    return ps;
}

void pauli_blocking_free(PauliState *ps) {
    if (ps) {
        free(ps->proj_targ_interaction);
        free(ps);
    }
}

void pauli_blocking_mark_interaction(PauliState *ps, int particle_idx) {
    if (ps && ps->proj_targ_interaction) {
        ps->proj_targ_interaction[particle_idx] = 1;
    }
}

int pauli_blocking_check(const SimState *state,
                          PauliState *ps,
                          const Particle *p1,
                          const Particle *p2) {
    if (!ps || !state || !p1 || !p2) return 0;

    const SimConfig *cfg = state->config;
    const Grid *grid = &cfg->grid;
    const Particle *particles = state->particles;
    float D_pot = cfg->potential.D;
    int num_ensembles = cfg->num_ensembles;
    int pauli_ensemble_count = cfg->collision.pauli_ensemble_count;

    /* Find particle indices */
    int idx1 = p1 - particles;
    int idx2 = p2 - particles;

    /* Cold nucleus filter: check if both particles are from same cold nucleus
     * If both from projectile and neither has hit target, block
     * If both from target and neither has hit projectile, block */
    int nucleon1_type = p1->is_projectile;
    int nucleon2_type = p2->is_projectile;

    if (nucleon1_type == nucleon2_type) {
        if (!ps->proj_targ_interaction[idx1] && !ps->proj_targ_interaction[idx2]) {
            return 1;  /* BLOCKED: cold nucleus collision */
        }
    }

    /* Check Pauli blocking for both particles in collision */
    for (int two_particle = 1; two_particle <= 2; two_particle++) {
        const Particle *p_local;
        int isospin_species;
        int current_ensemble;

        if (two_particle == 1) {
            p_local = p1;
            isospin_species = p1->charge;
            current_ensemble = p1->ensemble_id;
        } else {
            p_local = p2;
            isospin_species = p2->charge;
            current_ensemble = p2->ensemble_id;
        }

        /* Local Pauli parameters (may be adjusted for isospin) */
        float r_pauli2_local = ps->r_pauli2;
        float p_pauli2_local = ps->p_pauli2;
        float r_pauli_local = ps->r_pauli;
        float p_pauli_local = ps->p_pauli;
        float local_occupancy = ps->max_occupancy;

        /* Apply isospin-dependent radii */
        if (D_pot != 0.0) {
            if (isospin_species) {  /* proton */
                r_pauli2_local = ps->r_pauli2 * square(ps->r_proton_iso_factor);
                r_pauli_local = ps->r_pauli * ps->r_proton_iso_factor;
                p_pauli2_local = ps->p_pauli2 * square(ps->p_proton_iso_factor);
                p_pauli_local = ps->p_pauli * ps->p_proton_iso_factor;
            } else {  /* neutron */
                r_pauli2_local = ps->r_pauli2 * square(ps->r_neutron_iso_factor);
                r_pauli_local = ps->r_pauli * ps->r_neutron_iso_factor;
                p_pauli2_local = ps->p_pauli2 * square(ps->p_neutron_iso_factor);
                p_pauli_local = ps->p_pauli * ps->p_neutron_iso_factor;
            }
        }

        /* Select ensembles to sample */
        int *ensemble_selected = (int*)calloc(num_ensembles + 1, sizeof(int));
        if (!ensemble_selected) return 0;

        int n_pauli_ensembles = pauli_ensemble_count;
        if (n_pauli_ensembles > num_ensembles) {
            n_pauli_ensembles = num_ensembles;
        }

        if (n_pauli_ensembles == num_ensembles) {
            /* Use all ensembles */
            for (int j = 1; j <= n_pauli_ensembles; j++) {
                ensemble_selected[j] = 1;
            }
        } else {
            /* Random ensemble selection - always include current ensemble */
            ensemble_selected[current_ensemble] = 1;

            for (int j = 1; j < n_pauli_ensembles; j++) {
                int ensemble_select = (int)(random_number(1) * num_ensembles) + 1;
                if (ensemble_select > num_ensembles) {
                    ensemble_select = num_ensembles;
                }

                while (ensemble_selected[ensemble_select]) {
                    ensemble_select = (int)(random_number(1) * num_ensembles) + 1;
                    if (ensemble_select > num_ensembles) {
                        ensemble_select = num_ensembles;
                    }
                }
                ensemble_selected[ensemble_select] = 1;
            }
        }

        /* Count particles in same phase space cell */
        float blocking_count = 0.0;

        /* Find grid cell range to search */
        float dx_inv = 1.0 / grid->dx;
        int Lx_low = where_am_I(p_local->x - r_pauli_local, grid->origin_x, dx_inv);
        int Ly_low = where_am_I(p_local->y - r_pauli_local, grid->origin_y, dx_inv);
        int Lz_low = where_am_I(p_local->z - r_pauli_local, grid->origin_z, dx_inv);
        int Lx_hi = where_am_I(p_local->x + r_pauli_local, grid->origin_x, dx_inv);
        int Ly_hi = where_am_I(p_local->y + r_pauli_local, grid->origin_y, dx_inv);
        int Lz_hi = where_am_I(p_local->z + r_pauli_local, grid->origin_z, dx_inv);

        /* Loop over grid cells
         * NOTE: OpenMP not used here due to blocking_count race condition
         * and early-exit logic. The outer collision loop should be parallelized instead. */
        for (int i_grid = Lx_low; i_grid <= Lx_hi; i_grid++) {
            for (int j_grid = Ly_low; j_grid <= Ly_hi; j_grid++) {
                for (int k_grid = Lz_low; k_grid <= Lz_hi; k_grid++) {
                    /* Check valid grid site (1-indexed in legacy code) */
                    if (i_grid < 1 || i_grid > grid->nx ||
                        j_grid < 1 || j_grid > grid->ny ||
                        k_grid < 1 || k_grid > grid->nz) {
                        continue;
                    }

                    /* Convert to 0-indexed for modern grid */
                    int cell_idx = grid_index(grid, i_grid - 1, j_grid - 1, k_grid - 1);
                    int n_particles_in_cell = grid->cell_counts[cell_idx];

                    /* Loop over particles in this cell */
                    for (int n_grid = 0; n_grid < n_particles_in_cell; n_grid++) {
                        int id = grid->cell_particles[cell_idx][n_grid];
                        const Particle *p_test = &particles[id];

                        /* Check if this ensemble is selected */
                        if (!ensemble_selected[p_test->ensemble_id]) continue;

                        /* Check for isospin match (if D_pot != 0) */
                        if (D_pot != 0.0 && isospin_species != p_test->charge) {
                            continue;
                        }

                        /* Check momentum space separation (fast rejection) */
                        float delta_x = p_local->px - p_test->px;
                        float delta_y = p_local->py - p_test->py;
                        float delta_z = p_local->pz - p_test->pz;

                        if (fabsf(delta_x) > p_pauli_local) continue;
                        if (fabsf(delta_y) > p_pauli_local) continue;
                        if (fabsf(delta_z) > p_pauli_local) continue;

                        float momentum_check = square(delta_x) + square(delta_y) + square(delta_z);
                        if (momentum_check > p_pauli2_local) continue;

                        /* Check coordinate space separation */
                        delta_x = p_local->x - p_test->x;
                        delta_y = p_local->y - p_test->y;
                        delta_z = p_local->z - p_test->z;

                        if (fabsf(delta_x) > r_pauli_local) continue;
                        if (fabsf(delta_y) > r_pauli_local) continue;
                        if (fabsf(delta_z) > r_pauli_local) continue;

                        float radius_check = square(delta_x) + square(delta_y) + square(delta_z);
                        if (radius_check > r_pauli2_local) continue;

                        /* Particle is in same phase space element */
                        blocking_count++;

                        /* Early exit if definitely blocked */
                        if (blocking_count > (local_occupancy + 1.0)) {
                            free(ensemble_selected);
                            return 1;  /* BLOCKED */
                        }
                    }
                }
            }
        }

        /* Check phase space density (subtract 1 for the particle itself) */
        float blocking_check = (blocking_count - 1.0) / local_occupancy;
        if (blocking_check > random_number(1)) {
            free(ensemble_selected);
            return 1;  /* BLOCKED */
        }

        free(ensemble_selected);
    }

    return 0;  /* ALLOWED */
}
