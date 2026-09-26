#include "integrator.h"
#include "potential.h"
#include "grid.h"
#include "utils.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void integrator_calculate_forces(SimState *state) {
    /* Parallelize force calculation across particles */
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < state->num_particles; i++) {
        Particle *p = &state->particles[i];
        potential_calculate_force(&state->config->grid, &state->config->potential,
                                   p, &p->dpxdt, &p->dpydt, &p->dpzdt);
    }
}

void integrator_update_positions(SimState *state) {
    float dt = state->config->dt;

    /* Parallelize + vectorize position updates */
    #pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < state->num_particles; i++) {
        Particle *p = &state->particles[i];

        /* Relativistic velocity: v = p/E */
        float E = energy_relativistic(p->px, p->py, p->pz, NUCLEON_MASS);
        p->dxdt = p->px / E;
        p->dydt = p->py / E;
        p->dzdt = p->pz / E;

        /* Update positions */
        p->x += p->dxdt * dt;
        p->y += p->dydt * dt;
        p->z += p->dzdt * dt;
    }
}

void integrator_update_momenta(SimState *state) {
    float dt = state->config->dt;

    /* Parallelize + vectorize momentum updates */
    #pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < state->num_particles; i++) {
        Particle *p = &state->particles[i];

        /* Update momenta from forces */
        p->px += p->dpxdt * dt;
        p->py += p->dpydt * dt;
        p->pz += p->dpzdt * dt;
    }
}

void integrator_step(SimState *state) {
    /* Velocity Verlet integration */
    float dt = state->config->dt;

    /* Calculate forces at t */
    integrator_calculate_forces(state);

    /* Half-step momentum update */
    for (size_t i = 0; i < state->num_particles; i++) {
        Particle *p = &state->particles[i];
        p->px += 0.5f * p->dpxdt * dt;
        p->py += 0.5f * p->dpydt * dt;
        p->pz += 0.5f * p->dpzdt * dt;
    }

    /* Full-step position update */
    integrator_update_positions(state);

    /* Update grid densities and potentials */
    grid_assign_particles(&state->config->grid, state->particles, state->num_particles);
    grid_calculate_density(&state->config->grid, state->particles, state->num_particles,
                           state->config->form_factor_range);
    potential_calculate_meanfield(&state->config->grid, &state->config->potential);
    if (state->config->potential.enable_coulomb) {
        grid_solve_coulomb(&state->config->grid);
    }

    /* Calculate forces at t+dt */
    integrator_calculate_forces(state);

    /* Half-step momentum update */
    for (size_t i = 0; i < state->num_particles; i++) {
        Particle *p = &state->particles[i];
        p->px += 0.5f * p->dpxdt * dt;
        p->py += 0.5f * p->dpydt * dt;
        p->pz += 0.5f * p->dpzdt * dt;
    }

    state->current_timestep++;
    state->current_time += dt;
}
