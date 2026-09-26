#include "collision.h"
#include "utils.h"
#include <math.h>

/* Cugnon NN cross section parametrization (mb) */
static float sigma_nn_cugnon(float sqrt_s, int is_same_isospin) {
    /* Simplified version - constant for now */
    if (is_same_isospin) {
        return 40.0f;  /* pp or nn */
    } else {
        return 45.0f;  /* pn */
    }
}

int collision_check_geometric(const Particle *p1, const Particle *p2, float sigma_max) {
    /* Check if particles overlap cross-section disk */
    float dx = p2->x - p1->x;
    float dy = p2->y - p1->y;
    float dz = p2->z - p1->z;
    float r2 = dx*dx + dy*dy + dz*dz;

    /* Geometric cross section radius */
    float sigma_geom = sigma_max / PI;  /* mb -> fm^2, then radius */
    return r2 < sigma_geom;
}

float collision_cross_section(const Particle *p1, const Particle *p2, const CollisionParams *params) {
    /* Calculate CM energy */
    float E1 = energy_relativistic(p1->px, p1->py, p1->pz, NUCLEON_MASS);
    float E2 = energy_relativistic(p2->px, p2->py, p2->pz, NUCLEON_MASS);
    float s = (E1 + E2)*(E1 + E2) - (p1->px + p2->px)*(p1->px + p2->px)
                                    - (p1->py + p2->py)*(p1->py + p2->py)
                                    - (p1->pz + p2->pz)*(p1->pz + p2->pz);
    float sqrt_s = sqrtf(s);

    int same_isospin = (p1->charge == p2->charge);
    return sigma_nn_cugnon(sqrt_s, same_isospin);
}

int collision_scatter(Particle *p1, Particle *p2, const CollisionParams *params) {
    /* Get cross section */
    float sigma = collision_cross_section(p1, p2, params);

    /* Monte Carlo acceptance */
    if (random_uniform() > sigma / params->sigma_nn_max) {
        return 0;  /* Rejected */
    }

    /* Transform to CM frame */
    float px_cm = p1->px + p2->px;
    float py_cm = p1->py + p2->py;
    float pz_cm = p1->pz + p2->pz;
    float E_cm = energy_relativistic(p1->px, p1->py, p1->pz, NUCLEON_MASS) +
                 energy_relativistic(p2->px, p2->py, p2->pz, NUCLEON_MASS);

    float beta_x = px_cm / E_cm;
    float beta_y = py_cm / E_cm;
    float beta_z = pz_cm / E_cm;

    /* Boost to CM */
    float p1_px = p1->px, p1_py = p1->py, p1_pz = p1->pz;
    float p1_E = energy_relativistic(p1->px, p1->py, p1->pz, NUCLEON_MASS);
    lorentz_boost(&p1_px, &p1_py, &p1_pz, &p1_E, -beta_x, -beta_y, -beta_z);

    /* Isotropic scattering in CM */
    float p_mag = sqrtf(p1_px*p1_px + p1_py*p1_py + p1_pz*p1_pz);
    float cos_theta = 2.0f * random_uniform() - 1.0f;
    float sin_theta = sqrtf(1.0f - cos_theta*cos_theta);
    float phi = 2.0f * PI * random_uniform();

    float p1_new_px = p_mag * sin_theta * cosf(phi);
    float p1_new_py = p_mag * sin_theta * sinf(phi);
    float p1_new_pz = p_mag * cos_theta;

    /* Momentum conservation in CM */
    float p2_new_px = -p1_new_px;
    float p2_new_py = -p1_new_py;
    float p2_new_pz = -p1_new_pz;

    /* Boost back to lab */
    float p1_new_E = p1_E;
    float p2_new_E = energy_relativistic(p2->px, p2->py, p2->pz, NUCLEON_MASS);
    lorentz_boost(&p1_new_px, &p1_new_py, &p1_new_pz, &p1_new_E, beta_x, beta_y, beta_z);
    lorentz_boost(&p2_new_px, &p2_new_py, &p2_new_pz, &p2_new_E, beta_x, beta_y, beta_z);

    /* Update momenta */
    p1->px = p1_new_px;
    p1->py = p1_new_py;
    p1->pz = p1_new_pz;
    p2->px = p2_new_px;
    p2->py = p2_new_py;
    p2->pz = p2_new_pz;

    p1->num_collisions++;
    p2->num_collisions++;

    return 1;  /* Scattered */
}

int collision_check_pauli(const SimState *state, const Particle *p1, const Particle *p2) {
    /* Simplified Pauli blocking - always allow for now */
    /* Full implementation would check phase space occupancy */
    return 1;
}

int collision_process_all(SimState *state) {
    int num_collisions = 0;

    /* Loop over ensembles */
    for (int ens = 0; ens < state->config->num_ensembles; ens++) {
        /* Find particles in this ensemble */
        for (size_t i = 0; i < state->num_particles; i++) {
            if (state->particles[i].ensemble_id != ens) continue;

            for (size_t j = i+1; j < state->num_particles; j++) {
                if (state->particles[j].ensemble_id != ens) continue;

                Particle *p1 = &state->particles[i];
                Particle *p2 = &state->particles[j];

                /* Geometric test */
                if (!collision_check_geometric(p1, p2, state->config->collision.sigma_nn_max)) {
                    continue;
                }

                /* Scatter */
                if (collision_scatter(p1, p2, &state->config->collision)) {
                    /* Check Pauli blocking */
                    if (collision_check_pauli(state, p1, p2)) {
                        num_collisions++;
                    }
                }
            }
        }
    }

    state->total_collisions += num_collisions;
    return num_collisions;
}
