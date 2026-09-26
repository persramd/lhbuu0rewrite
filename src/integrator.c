#include "integrator.h"
#include "potential.h"
#include "grid.h"
#include "utils.h"
#include "constants.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Temporary storage for modified Velocity Verlet algorithm */
static float *p_half_x = NULL;
static float *p_half_y = NULL;
static float *p_half_z = NULL;
static float *a_zero_x = NULL;
static float *a_zero_y = NULL;
static float *a_zero_z = NULL;
static size_t storage_capacity = 0;

/* Allocate temporary storage for Verlet integration */
static void ensure_storage(size_t n_particles) {
    if (n_particles > storage_capacity) {
        free(p_half_x); free(p_half_y); free(p_half_z);
        free(a_zero_x); free(a_zero_y); free(a_zero_z);

        p_half_x = malloc(n_particles * sizeof(float));
        p_half_y = malloc(n_particles * sizeof(float));
        p_half_z = malloc(n_particles * sizeof(float));
        a_zero_x = malloc(n_particles * sizeof(float));
        a_zero_y = malloc(n_particles * sizeof(float));
        a_zero_z = malloc(n_particles * sizeof(float));

        storage_capacity = n_particles;
    }
}

/* Calculate forces: position-space gradient of potentials */
void integrator_calculate_forces(SimState *state) {
    const Grid *grid = &state->config->grid;
    const PotentialParams *pot = &state->config->potential;
    const int ff_range = state->config->form_factor_range;

    /* Parallelize force calculation across particles */
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < state->num_particles; i++) {
        Particle *p = &state->particles[i];

        /* Initialize forces to zero */
        float a_part[3] = {0.0f, 0.0f, 0.0f};

        /* Get particle grid location */
        int Lx, Ly, Lz;
        grid_coords(grid, p->x, p->y, p->z, &Lx, &Ly, &Lz);

        /* Check bounds (need ff_range buffer) */
        int min_grid = ff_range + 1;
        int max_x_grid = grid->nx - ff_range;
        int max_y_grid = grid->ny - ff_range;
        int max_z_grid = grid->nz - ff_range;

        if (Lx > min_grid && Lx < max_x_grid &&
            Ly > min_grid && Ly < max_y_grid &&
            Lz > min_grid && Lz < max_z_grid) {

            /* Calculate force from mean-field potentials
             * This is half_impulse() from legacy: calculates -∇U
             * using form factor gradients drho_dr */
            potential_calculate_force(grid, pot, p, &a_part[0], &a_part[1], &a_part[2]);
        }

        /* Store forces (accelerations) in particle structure */
        p->dpxdt = a_part[0];
        p->dpydt = a_part[1];
        p->dpzdt = a_part[2];
    }
}

/* Calculate momentum-dependent potential gradients for position updates */
static void calculate_momentum_gradients(SimState *state, float *p_grad_x, float *p_grad_y, float *p_grad_z) {
    const Grid *grid = &state->config->grid;
    const PotentialParams *pot = &state->config->potential;
    const int ff_range = state->config->form_factor_range;

    /* Initialize all gradients to zero */
    memset(p_grad_x, 0, state->num_particles * sizeof(float));
    memset(p_grad_y, 0, state->num_particles * sizeof(float));
    memset(p_grad_z, 0, state->num_particles * sizeof(float));

    /* Only compute if momentum-dependent potential is enabled */
    if (!pot->enable_momentum_dep || pot->C == 0.0f) {
        return;
    }

    /* Parallelize gradient calculation */
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < state->num_particles; i++) {
        Particle *p = &state->particles[i];

        /* Get particle grid location */
        int Lx, Ly, Lz;
        grid_coords(grid, p->x, p->y, p->z, &Lx, &Ly, &Lz);

        /* Check bounds */
        int min_grid = ff_range + 1;
        int max_x_grid = grid->nx - ff_range;
        int max_y_grid = grid->ny - ff_range;
        int max_z_grid = grid->nz - ff_range;

        if (Lx > min_grid && Lx < max_x_grid &&
            Ly > min_grid && Ly < max_y_grid &&
            Lz > min_grid && Lz < max_z_grid) {

            /* Calculate ∇_p V(r,p) for momentum-dependent potential
             * This is grad_V_alpha_p_md() from legacy */
            float gx, gy, gz;
            potential_grad_momentum_dep(grid, pot, p, Lx, Ly, Lz, &gx, &gy, &gz);

            p_grad_x[i] = gx;
            p_grad_y[i] = gy;
            p_grad_z[i] = gz;
        }
    }
}

/* Modified Velocity Verlet integration step
 * From legacy: int_the_great_eom.c (Allen & Tildesley algorithm)
 *
 * Key modification: position update includes momentum-dependent gradient
 * r(t+dt) = r(t) + p_half*dt/E + ∇_p V(r,p)*dt
 */
void integrator_step(SimState *state) {
    const float dt = state->config->dt;
    const size_t n = state->num_particles;

    ensure_storage(n);

    /* ===== STEP 1: Calculate forces at t ===== */
    integrator_calculate_forces(state);

    /* ===== STEP 2: Half-step momentum predictor =====
     * p_half = p(t) + 0.5*F(t)*dt
     * Store both p_half and full-step predictor a_zero = 2*F(t)*dt/2 = F(t)*dt */
    #pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < n; i++) {
        Particle *p = &state->particles[i];

        /* Half-step acceleration */
        float a_x = p->dpxdt * 0.5f * dt;
        float a_y = p->dpydt * 0.5f * dt;
        float a_z = p->dpzdt * 0.5f * dt;

        /* Store half-step momentum for position update */
        p_half_x[i] = p->px + a_x;
        p_half_y[i] = p->py + a_y;
        p_half_z[i] = p->pz + a_z;

        /* Store full-step accelerations (used later for momentum-dependent potential) */
        a_zero_x[i] = 2.0f * a_x;
        a_zero_y[i] = 2.0f * a_y;
        a_zero_z[i] = 2.0f * a_z;
    }

    /* ===== STEP 3: Calculate momentum-dependent gradients at p_half ===== */
    float *p_grad_x = malloc(n * sizeof(float));
    float *p_grad_y = malloc(n * sizeof(float));
    float *p_grad_z = malloc(n * sizeof(float));

    calculate_momentum_gradients(state, p_grad_x, p_grad_y, p_grad_z);

    /* ===== STEP 4: Full-step position update with momentum-dependent term =====
     * r(t+dt) = r(t) + v*dt + ∇_p V(r,p)*dt
     * where v = p_half/E (relativistic velocity) */
    #pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < n; i++) {
        Particle *p = &state->particles[i];

        /* Relativistic energy at half-step momentum */
        float E = sqrtf(p_half_x[i]*p_half_x[i] +
                       p_half_y[i]*p_half_y[i] +
                       p_half_z[i]*p_half_z[i] +
                       NUCLEON_MASS*NUCLEON_MASS);
        float dt_over_E = dt / E;

        /* Update positions: r += (p_half/E)*dt + ∇_p V*dt */
        p->x += p_half_x[i] * dt_over_E + p_grad_x[i] * dt;
        p->y += p_half_y[i] * dt_over_E + p_grad_y[i] * dt;
        p->z += p_half_z[i] * dt_over_E + p_grad_z[i] * dt;

        /* Store velocities for diagnostics */
        p->dxdt = p_half_x[i] / E + p_grad_x[i];
        p->dydt = p_half_y[i] / E + p_grad_y[i];
        p->dzdt = p_half_z[i] / E + p_grad_z[i];
    }

    free(p_grad_x); free(p_grad_y); free(p_grad_z);

    /* ===== STEP 5: Update grid densities and potentials at new positions ===== */
    grid_assign_particles(&state->config->grid, state->particles, n);
    grid_calculate_density(&state->config->grid, state->particles, n,
                          state->config->form_factor_range);
    potential_calculate_meanfield(&state->config->grid, &state->config->potential);

    if (state->config->potential.enable_coulomb) {
        grid_solve_coulomb(&state->config->grid);
    }

    /* ===== STEP 6: For momentum-dependent potential, approximate p(t+dt) ===== */
    if (state->config->potential.enable_momentum_dep && state->config->potential.C != 0.0f) {
        /* Temporarily advance momenta to full-step for MDYI evaluation
         * This is needed because U(r,p) depends on current momentum */
        #pragma omp parallel for simd schedule(static)
        for (size_t i = 0; i < n; i++) {
            Particle *p = &state->particles[i];
            p->px += a_zero_x[i];
            p->py += a_zero_y[i];
            p->pz += a_zero_z[i];
        }

        /* TODO: Update momentum-dependent potential grid here (p_avg_alpha) */
    }

    /* ===== STEP 7: Calculate forces at t+dt ===== */
    integrator_calculate_forces(state);

    /* ===== STEP 8: Half-step momentum corrector =====
     * p(t+dt) = p_half + 0.5*F(t+dt)*dt */
    #pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < n; i++) {
        Particle *p = &state->particles[i];

        /* Final momentum: p_half + second half-step */
        p->px = p_half_x[i] + p->dpxdt * 0.5f * dt;
        p->py = p_half_y[i] + p->dpydt * 0.5f * dt;
        p->pz = p_half_z[i] + p->dpzdt * 0.5f * dt;
    }

    /* TODO: Update momentum-dependent potential grid with final momenta (p_avg_alpha) */

    /* Advance time */
    state->current_timestep++;
    state->current_time += dt;
}

/* Cleanup storage on exit */
void integrator_cleanup(void) {
    free(p_half_x); free(p_half_y); free(p_half_z);
    free(a_zero_x); free(a_zero_y); free(a_zero_z);
    p_half_x = p_half_y = p_half_z = NULL;
    a_zero_x = a_zero_y = a_zero_z = NULL;
    storage_capacity = 0;
}
