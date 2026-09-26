#include "collision.h"
#include "grid.h"
#include "utils.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Cugnon NN cross section (NIM B111, 1996) - full parametrization */
static float sigma_nn_cugnon_full(float sqrt_s, int is_pp, int is_nn, int is_pn) {
    /* Convert to lab kinetic energy */
    float m = NUCLEON_MASS;
    float T_lab = sqrt_s*sqrt_s / (2.0f * m) - m;  /* GeV */
    T_lab *= 1000.0f;  /* MeV */

    float sigma;

    if (T_lab < 0.4f) {
        /* Low energy - use constant */
        if (is_pp || is_nn) {
            sigma = 10.0f;  /* mb */
        } else {
            sigma = 33.0f;  /* pn */
        }
    } else if (T_lab < 2000.0f) {
        /* Cugnon parametrization */
        float a, b, c, d, e, f;

        if (is_pp || is_nn) {
            /* pp or nn */
            a = 35.80f;
            b = 0.308f;
            c = 17.0f;
            d = 0.0115f;
            e = 0.01f;
            f = 4.0f;
        } else {
            /* pn */
            a = 35.45f;
            b = 0.3f;
            c = 20.0f;
            d = 0.008f;
            e = 0.01f;
            f = 5.5f;
        }

        sigma = a + b * powf(T_lab, -c) + d * T_lab * T_lab - e * powf(T_lab, f);
    } else {
        /* High energy asymptotic */
        sigma = 30.0f + 0.15f * logf(T_lab / 200.0f);
    }

    return fmaxf(sigma, 1.0f);  /* Minimum 1 mb */
}

/* Density-dependent cross section reduction */
static float sigma_medium_factor(float rho, float rho0, float reduction_factor) {
    return 1.0f - reduction_factor * (rho / rho0);
}

/* Relative velocity in NN CM frame */
static float relative_velocity(const Particle *p1, const Particle *p2) {
    /* Total momentum */
    float Px = p1->px + p2->px;
    float Py = p1->py + p2->py;
    float Pz = p1->pz + p2->pz;

    /* Total energy */
    float E1 = energy_relativistic(p1->px, p1->py, p1->pz, NUCLEON_MASS);
    float E2 = energy_relativistic(p2->px, p2->py, p2->pz, NUCLEON_MASS);
    float E_tot = E1 + E2;

    /* CM velocity */
    float beta_x = Px / E_tot;
    float beta_y = Py / E_tot;
    float beta_z = Pz / E_tot;

    /* Boost p1 to CM */
    float p1_cm_x = p1->px;
    float p1_cm_y = p1->py;
    float p1_cm_z = p1->pz;
    float E1_cm = E1;
    lorentz_boost(&p1_cm_x, &p1_cm_y, &p1_cm_z, &E1_cm,
                  -beta_x, -beta_y, -beta_z);

    /* Relative velocity in CM */
    float p_cm = sqrtf(p1_cm_x*p1_cm_x + p1_cm_y*p1_cm_y + p1_cm_z*p1_cm_z);
    float v_rel = 2.0f * p_cm / E1_cm;  /* Relativistic relative velocity */

    return v_rel;
}

/* Lorentz-contracted distance for collision check */
static float contracted_distance(const Particle *p1, const Particle *p2,
                                 float dx, float dy, float dz) {
    /* CM velocity */
    float Px = p1->px + p2->px;
    float Py = p1->py + p2->py;
    float Pz = p1->pz + p2->pz;

    float E1 = energy_relativistic(p1->px, p1->py, p1->pz, NUCLEON_MASS);
    float E2 = energy_relativistic(p2->px, p2->py, p2->pz, NUCLEON_MASS);
    float E_tot = E1 + E2;

    float beta_x = Px / E_tot;
    float beta_y = Py / E_tot;
    float beta_z = Pz / E_tot;
    float beta = sqrtf(beta_x*beta_x + beta_y*beta_y + beta_z*beta_z);

    if (beta < 1e-6f) {
        return sqrtf(dx*dx + dy*dy + dz*dz);
    }

    /* Parallel and perpendicular components */
    float beta_dot_r = (beta_x * dx + beta_y * dy + beta_z * dz) / beta;
    float r_perp_sq = dx*dx + dy*dy + dz*dz - beta_dot_r * beta_dot_r;

    /* Lorentz contraction */
    float gamma = 1.0f / sqrtf(1.0f - beta*beta);
    float r_parallel = beta_dot_r / gamma;

    return sqrtf(r_parallel*r_parallel + r_perp_sq);
}

/* Check Pauli blocking with proper phase-space occupancy */
static int check_pauli_blocking(const SimState *state,
                                const Particle *p1_new, const Particle *p2_new) {
    const Grid *grid = &state->config->grid;
    float pf = 0.27f;  /* Fermi momentum at rho0 (GeV) */

    /* Check both final-state particles */
    const Particle *particles[2] = {p1_new, p2_new};

    for (int i = 0; i < 2; i++) {
        const Particle *p = particles[i];

        /* Get grid cell */
        int ix, iy, iz;
        grid_coords(grid, p->x, p->y, p->z, &ix, &iy, &iz);

        if (ix < 0 || ix >= grid->nx || iy < 0 || iy >= grid->ny ||
            iz < 0 || iz >= grid->nz) {
            continue;  /* Outside grid - allow */
        }

        /* Get local density */
        int idx = grid_index(grid, ix, iy, iz);
        float rho = grid->baryon_density[idx];

        /* Local Fermi momentum */
        float pf_local = pf * powf(rho / 0.16f, 1.0f/3.0f);

        /* Check if momentum is above Fermi surface */
        float p_mag = sqrtf(p->px*p->px + p->py*p->py + p->pz*p->pz);

        if (p_mag < pf_local) {
            /* Below Fermi surface - check occupancy */
            /* Count particles in nearby phase space */
            int count_same_type = 0;
            int cell_idx = grid_index(grid, ix, iy, iz);

            for (int k = 0; k < grid->cell_counts[cell_idx]; k++) {
                int pidx = grid->cell_particles[cell_idx][k];
                const Particle *other = &state->particles[pidx];

                /* Same charge */
                if (other->charge != p->charge) continue;

                /* Similar momentum */
                float dp = sqrtf((other->px - p->px)*(other->px - p->px) +
                               (other->py - p->py)*(other->py - p->py) +
                               (other->pz - p->pz)*(other->pz - p->pz));

                if (dp < 0.1f) {  /* Momentum cell size */
                    count_same_type++;
                }
            }

            /* Pauli blocking if phase space occupied */
            if (count_same_type >= 4) {  /* Spin-isospin degeneracy */
                return 1;  /* Blocked */
            }
        }
    }

    return 0;  /* Allowed */
}

int collision_process_all(SimState *state) {
    int num_collisions = 0;
    int num_pauli_blocked = 0;

    /* Loop over ensembles */
    for (int ens = 0; ens < state->config->num_ensembles; ens++) {
        /* Find particles in this ensemble */
        for (size_t i = 0; i < state->num_particles; i++) {
            if (state->particles[i].ensemble_id != ens) continue;

            for (size_t j = i+1; j < state->num_particles; j++) {
                if (state->particles[j].ensemble_id != ens) continue;

                Particle *p1 = &state->particles[i];
                Particle *p2 = &state->particles[j];

                /* Spatial separation */
                float dx = p2->x - p1->x;
                float dy = p2->y - p1->y;
                float dz = p2->z - p1->z;

                /* Lorentz-contracted distance */
                float r_contract = contracted_distance(p1, p2, dx, dy, dz);

                /* Geometric test with Lorentz contraction */
                float sigma_geom = state->config->collision.sigma_nn_max / PI;
                if (r_contract > sqrtf(sigma_geom)) {
                    continue;  /* Too far */
                }

                /* Get collision properties */
                int is_pp = (p1->charge == 1 && p2->charge == 1);
                int is_nn = (p1->charge == 0 && p2->charge == 0);
                int is_pn = !is_pp && !is_nn;

                /* CM energy */
                float E1 = energy_relativistic(p1->px, p1->py, p1->pz, NUCLEON_MASS);
                float E2 = energy_relativistic(p2->px, p2->py, p2->pz, NUCLEON_MASS);
                float s = (E1 + E2)*(E1 + E2) -
                         ((p1->px + p2->px)*(p1->px + p2->px) +
                          (p1->py + p2->py)*(p1->py + p2->py) +
                          (p1->pz + p2->pz)*(p1->pz + p2->pz));
                float sqrt_s = sqrtf(fmaxf(s, 0.0f));

                /* Free NN cross section */
                float sigma_free = sigma_nn_cugnon_full(sqrt_s, is_pp, is_nn, is_pn);

                /* In-medium modification */
                int ix, iy, iz;
                grid_coords(&state->config->grid,
                           0.5f*(p1->x + p2->x),
                           0.5f*(p1->y + p2->y),
                           0.5f*(p1->z + p2->z),
                           &ix, &iy, &iz);

                float rho_local = 0.0f;
                if (ix >= 0 && ix < state->config->grid.nx &&
                    iy >= 0 && iy < state->config->grid.ny &&
                    iz >= 0 && iz < state->config->grid.nz) {
                    int idx = grid_index(&state->config->grid, ix, iy, iz);
                    rho_local = state->config->grid.baryon_density[idx];
                }

                float medium_factor = sigma_medium_factor(rho_local, 0.16f,
                                    state->config->collision.sigma_rho_reduction);
                float sigma_eff = sigma_free * medium_factor;

                /* Monte Carlo test */
                float prob = sigma_eff / state->config->collision.sigma_nn_max;
                if (random_uniform() > prob) {
                    continue;  /* Rejected */
                }

                /* Scatter - store original momenta */
                float p1_orig_px = p1->px, p1_orig_py = p1->py, p1_orig_pz = p1->pz;
                float p2_orig_px = p2->px, p2_orig_py = p2->py, p2_orig_pz = p2->pz;

                if (collision_scatter(p1, p2, &state->config->collision)) {
                    /* Check Pauli blocking */
                    if (state->config->collision.enable_pauli_blocking &&
                        check_pauli_blocking(state, p1, p2)) {
                        /* Blocked - restore original momenta */
                        p1->px = p1_orig_px;
                        p1->py = p1_orig_py;
                        p1->pz = p1_orig_pz;
                        p2->px = p2_orig_px;
                        p2->py = p2_orig_py;
                        p2->pz = p2_orig_pz;
                        num_pauli_blocked++;
                    } else {
                        /* Collision successful */
                        num_collisions++;
                    }
                }
            }
        }
    }

    state->total_collisions += num_collisions;
    state->pauli_blocked += num_pauli_blocked;

    return num_collisions;
}
