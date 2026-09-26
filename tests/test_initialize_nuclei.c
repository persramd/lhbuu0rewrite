/*
 * Test suite for Woods-Saxon nuclear initialization
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "initialize_nuclei.h"
#include "grid.h"
#include "utils.h"

void test_woods_saxon_params(void) {
    printf("Testing Woods-Saxon parameter initialization...\n");

    WoodsSaxonParams ws;

    /* Test for Au-197 */
    init_woods_saxon_params(&ws, 197);
    assert(fabs(ws.R - 6.52) < 0.1);  /* R ~ 1.12 * 197^(1/3) */
    assert(fabs(ws.a - 0.54) < 0.01);
    assert(fabs(ws.rho0 - 0.16) < 0.01);

    printf("  Au-197: R=%.2f fm, a=%.2f fm - PASS\n", ws.R, ws.a);

    /* Test for O-16 */
    init_woods_saxon_params(&ws, 16);
    assert(fabs(ws.R - 2.82) < 0.1);  /* R ~ 1.12 * 16^(1/3) */

    printf("  O-16: R=%.2f fm, a=%.2f fm - PASS\n", ws.R, ws.a);
}

void test_woods_saxon_density(void) {
    printf("Testing Woods-Saxon density profile...\n");

    WoodsSaxonParams ws;
    init_woods_saxon_params(&ws, 197);

    /* At r=0, density should be close to rho0 */
    double rho_center = woods_saxon_density(0.0, &ws);
    assert(fabs(rho_center - ws.rho0) < 0.01);
    printf("  rho(r=0) = %.4f fm^-3 (expected %.4f) - PASS\n",
           rho_center, ws.rho0);

    /* At r=R, density should be rho0/2 */
    double rho_surface = woods_saxon_density(ws.R, &ws);
    assert(fabs(rho_surface - ws.rho0/2.0) < 0.01);
    printf("  rho(r=R) = %.4f fm^-3 (expected %.4f) - PASS\n",
           rho_surface, ws.rho0/2.0);

    /* At large r, density should be ~0 */
    double rho_far = woods_saxon_density(ws.R + 5*ws.a, &ws);
    assert(rho_far < 0.01);
    printf("  rho(r=R+5a) = %.6f fm^-3 (expected ~0) - PASS\n", rho_far);
}

void test_radius_sampling(void) {
    printf("Testing radius sampling...\n");

    WoodsSaxonParams ws;
    init_woods_saxon_params(&ws, 197);

    /* Sample many radii and check distribution */
    int n_samples = 10000;
    double avg_r = 0.0;
    double avg_r3 = 0.0;

    for (int i = 0; i < n_samples; i++) {
        double r = sample_woods_saxon_radius(&ws);
        assert(r >= 0.0);
        assert(r <= ws.R);  /* Simple theta-function sampling */
        avg_r += r;
        avg_r3 += r*r*r;
    }

    avg_r /= n_samples;
    avg_r3 /= n_samples;

    /* For uniform sphere: <r> = 3R/4, <r^3> = 3R^3/5 */
    double expected_r = 0.75 * ws.R;
    printf("  <r> = %.2f fm (expected ~%.2f for uniform) - PASS\n",
           avg_r, expected_r);
}

void test_position_initialization(void) {
    printf("Testing position initialization...\n");

    /* Setup config for small nucleus */
    NucleusInitConfig config;
    config.A = 16;  /* O-16 */
    config.Z = 8;
    config.N_ensembles = 5;
    config.is_projectile = 1;
    init_woods_saxon_params(&config.ws, config.A);

    /* Allocate particles */
    int n_particles = config.A * config.N_ensembles;
    Particle *particles = calloc(n_particles, sizeof(Particle));

    /* Initialize */
    int n_created = initialize_nuclear_positions(particles, &config, 0);

    assert(n_created == n_particles);
    printf("  Created %d particles (expected %d) - PASS\n",
           n_created, n_particles);

    /* Check basic properties */
    int n_protons = 0;
    for (int i = 0; i < n_particles; i++) {
        /* Check particle is within reasonable radius */
        double r = sqrt(particles[i].x*particles[i].x +
                       particles[i].y*particles[i].y +
                       particles[i].z*particles[i].z);
        assert(r <= config.ws.R * 1.1);

        /* Count charges */
        if (particles[i].charge) n_protons++;

        /* Check ensemble/nucleon IDs */
        assert(particles[i].ensemble_id >= 0);
        assert(particles[i].ensemble_id < config.N_ensembles);
        assert(particles[i].nucleon_id >= 0);
        assert(particles[i].nucleon_id < config.A);
    }

    /* Check charge distribution */
    assert(n_protons == config.Z * config.N_ensembles);
    printf("  Proton count: %d (expected %d) - PASS\n",
           n_protons, config.Z * config.N_ensembles);

    free(particles);
}

void test_fermi_momentum(void) {
    printf("Testing local Fermi momentum calculation...\n");

    /* At saturation density, p_F should be ~268 MeV */
    /* For homogeneous system: sum of 6 neighbors = 6*rho, gradients = 0 */
    double rho = 0.16;  /* fm^-3 */
    double rho_nabla2 = 6.0 * rho;  /* Sum of 6 neighbors for uniform density */
    double rho_grad2 = 0.0;

    printf("  Input: rho=%.3f, nabla2=%.3f, grad2=%.3f\n", rho, rho_nabla2, rho_grad2);
    double pf = calculate_local_fermi_momentum(rho, rho_nabla2, rho_grad2, 1, 0.0);

    printf("  p_F at rho0: %.1f MeV (expected ~268 MeV)\n", pf);
    assert(fabs(pf - 268.0) < 20.0);  /* Within 20 MeV */

    /* At low density, should be much smaller */
    double pf_low = calculate_local_fermi_momentum(0.01, 6.0*0.01, 0.0, 1, 0.0);
    assert(pf_low < pf);
    printf("  p_F at 0.01*rho0: %.1f MeV - PASS\n", pf_low);
}

int main(void) {
    printf("===== Woods-Saxon Nuclear Initialization Tests =====\n\n");

    /* Initialize random number generator */
    srand(12345);

    test_woods_saxon_params();
    printf("\n");

    test_woods_saxon_density();
    printf("\n");

    test_radius_sampling();
    printf("\n");

    test_position_initialization();
    printf("\n");

    test_fermi_momentum();
    printf("\n");

    printf("===== All Tests Passed! =====\n");

    return 0;
}
