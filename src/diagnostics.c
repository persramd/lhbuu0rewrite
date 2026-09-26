#include "diagnostics.h"
#include "grid.h"
#include "utils.h"
#include <math.h>

float diagnostics_kinetic_energy(const SimState *state) {
    float ke = 0.0f;
    for (size_t i = 0; i < state->num_particles; i++) {
        const Particle *p = &state->particles[i];
        float E = energy_relativistic(p->px, p->py, p->pz, NUCLEON_MASS);
        ke += (E - NUCLEON_MASS);  /* Kinetic only */
    }
    return ke;
}

float diagnostics_potential_energy(const SimState *state) {
    float pe = 0.0f;
    const Grid *grid = &state->config->grid;

    for (size_t i = 0; i < state->num_particles; i++) {
        const Particle *p = &state->particles[i];
        int ix, iy, iz;
        grid_coords(grid, p->x, p->y, p->z, &ix, &iy, &iz);

        if (ix >= 0 && ix < grid->nx && iy >= 0 && iy < grid->ny && iz >= 0 && iz < grid->nz) {
            int idx = grid_index(grid, ix, iy, iz);
            pe += grid->potential[idx];
            if (state->config->potential.enable_coulomb && p->charge) {
                pe += grid->coulomb_potential[idx];
            }
        }
    }
    return pe / 2.0f;  /* Avoid double counting */
}

float diagnostics_total_energy(const SimState *state) {
    return diagnostics_kinetic_energy(state) + diagnostics_potential_energy(state);
}

void diagnostics_momentum(const SimState *state, float *px, float *py, float *pz) {
    *px = *py = *pz = 0.0f;
    for (size_t i = 0; i < state->num_particles; i++) {
        const Particle *p = &state->particles[i];
        *px += p->px;
        *py += p->py;
        *pz += p->pz;
    }
}

void diagnostics_angular_momentum(const SimState *state, float *lx, float *ly, float *lz) {
    *lx = *ly = *lz = 0.0f;
    for (size_t i = 0; i < state->num_particles; i++) {
        const Particle *p = &state->particles[i];
        *lx += p->y * p->pz - p->z * p->py;
        *ly += p->z * p->px - p->x * p->pz;
        *lz += p->x * p->py - p->y * p->px;
    }
}

void diagnostics_write_timestep(FILE *fp, const SimState *state) {
    float ke = diagnostics_kinetic_energy(state);
    float pe = diagnostics_potential_energy(state);
    float px, py, pz;
    diagnostics_momentum(state, &px, &py, &pz);

    fprintf(fp, "%d %.3f %.6f %.6f %.6f %.6f %.6f %.6f %d\n",
            state->current_timestep,
            state->current_time,
            ke, pe, ke + pe,
            px, py, pz,
            state->total_collisions);
}

int diagnostics_write_phase_space(const char *filename, const SimState *state) {
    FILE *fp = file_open(filename, "w");

    fprintf(fp, "# LHBUU Phase Space Output\n");
    fprintf(fp, "# Time: %.3f fm/c\n", state->current_time);
    fprintf(fp, "# Particles: %zu\n", state->num_particles);
    fprintf(fp, "# x y z px py pz charge ensemble\n");

    for (size_t i = 0; i < state->num_particles; i++) {
        const Particle *p = &state->particles[i];
        fprintf(fp, "%.4f %.4f %.4f %.6f %.6f %.6f %d %d\n",
                p->x, p->y, p->z,
                p->px, p->py, p->pz,
                p->charge, p->ensemble_id);
    }

    fclose(fp);
    return 1;
}

void diagnostics_print_summary(const SimState *state) {
    printf("\n=== LHBUU Simulation Summary ===\n");
    printf("Final time: %.2f fm/c\n", state->current_time);
    printf("Timesteps: %d\n", state->current_timestep);
    printf("Total collisions: %d\n", state->total_collisions);
    printf("Kinetic energy: %.4f GeV\n", diagnostics_kinetic_energy(state));
    printf("Potential energy: %.4f GeV\n", diagnostics_potential_energy(state));
    printf("Total energy: %.4f GeV\n", diagnostics_total_energy(state));

    float px, py, pz;
    diagnostics_momentum(state, &px, &py, &pz);
    printf("Momentum: (%.6f, %.6f, %.6f) GeV\n", px, py, pz);
    printf("================================\n\n");
}
