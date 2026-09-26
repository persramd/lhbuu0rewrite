#include "initialize.h"
#include "grid.h"
#include "potential.h"
#include "utils.h"
#include <stdlib.h>
#include <math.h>

void initialize_woods_saxon(Nucleus *nucleus) {
    /* Woods-Saxon parameters */
    float r0 = 1.12f * powf(nucleus->A, 1.0f/3.0f);  /* fm */
    float a = 0.54f;  /* fm, surface diffuseness */

    nucleus->radius = r0;

    /* For simplicity, store constant profile */
    nucleus->density_profile = NULL;
}

int initialize_nuclei(SimState *state) {
    SimConfig *cfg = state->config;

    /* Calculate nuclear radii */
    initialize_woods_saxon(&cfg->projectile);
    initialize_woods_saxon(&cfg->target);

    /* Calculate number of test particles */
    int n_proj = cfg->projectile.A * cfg->num_ensembles;
    int n_targ = cfg->target.A * cfg->num_ensembles;
    state->num_particles = n_proj + n_targ;
    state->max_particles = state->num_particles;

    state->particles = calloc(state->num_particles, sizeof(Particle));
    if (!state->particles) {
        error_exit("Failed to allocate particles");
    }

    /* Generate projectile particles */
    size_t idx = 0;
    for (int ens = 0; ens < cfg->num_ensembles; ens++) {
        for (int i = 0; i < cfg->projectile.A; i++) {
            Particle *p = &state->particles[idx++];

            /* Random position in sphere */
            float r, cos_theta, phi;
            do {
                r = cfg->projectile.radius * powf(random_uniform(), 1.0f/3.0f);
                cos_theta = 2.0f * random_uniform() - 1.0f;
                phi = 2.0f * PI * random_uniform();
            } while (0);

            float sin_theta = sqrtf(1.0f - cos_theta*cos_theta);
            p->x = r * sin_theta * cosf(phi);
            p->y = r * sin_theta * sinf(phi);
            p->z = r * cos_theta - 10.0f;  /* Offset projectile */

            p->px = 0.0f;
            p->py = 0.0f;
            p->pz = sqrtf(2.0f * NUCLEON_MASS * cfg->beam_energy / 1000.0f);  /* MeV->GeV */

            p->charge = (i < cfg->projectile.Z) ? 1 : 0;
            p->ensemble_id = ens;
            p->nucleon_id = i;
            p->is_projectile = 1;
        }
    }

    /* Generate target particles */
    for (int ens = 0; ens < cfg->num_ensembles; ens++) {
        for (int i = 0; i < cfg->target.A; i++) {
            Particle *p = &state->particles[idx++];

            /* Random position in sphere */
            float r = cfg->target.radius * powf(random_uniform(), 1.0f/3.0f);
            float cos_theta = 2.0f * random_uniform() - 1.0f;
            float phi = 2.0f * PI * random_uniform();

            float sin_theta = sqrtf(1.0f - cos_theta*cos_theta);
            float b = cfg->impact_parameter * (cfg->projectile.radius + cfg->target.radius);
            p->x = r * sin_theta * cosf(phi) + b;
            p->y = r * sin_theta * sinf(phi);
            p->z = r * cos_theta;

            p->px = 0.0f;
            p->py = 0.0f;
            p->pz = 0.0f;

            p->charge = (i < cfg->target.Z) ? 1 : 0;
            p->ensemble_id = ens;
            p->nucleon_id = i;
            p->is_projectile = 0;
        }
    }

    return 1;
}

int sim_init(SimState *state, SimConfig *config) {
    state->config = config;
    state->current_timestep = 0;
    state->current_time = 0.0f;
    state->total_collisions = 0;
    state->pauli_blocked = 0;

    /* Initialize grid */
    if (!grid_init(&config->grid, config->grid.nx, config->grid.ny, config->grid.nz, config->grid.dx)) {
        return 0;
    }

    /* Initialize particles */
    if (!initialize_nuclei(state)) {
        return 0;
    }

    /* Initial grid setup */
    grid_assign_particles(&config->grid, state->particles, state->num_particles);
    grid_calculate_density(&config->grid, state->particles, state->num_particles, config->form_factor_range);
    potential_calculate_meanfield(&config->grid, &config->potential);

    return 1;
}

void sim_free(SimState *state) {
    free(state->particles);
    grid_free(&state->config->grid);
}
