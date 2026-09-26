/*
 * LHBUU Lorentz Transformation Utilities - Implementation
 *
 * Ported from original code:
 *   /Users/declan/Documents/dev/osx/lhbuu/lib/lorentz_boostnntoNN.c
 *   /Users/declan/Documents/dev/osx/lhbuu/lib/collisions.c (lorentz() function)
 *
 * Key formulas:
 *   CM velocity: beta_i = (p1_i + p2_i) / (E1 + E2)
 *   CM energy: E_cm = sqrt((E1+E2)^2 - |p1+p2|^2)
 *   Boost: p'_parallel = (gamma-1)/beta^2 * (p·beta) - gamma*E
 *          p'_vector = p + p'_parallel * beta
 */

#include "lorentz.h"
#include "macros.h"
#include "types.h"
#include <math.h>
#include <stddef.h>

/* ================================================================
 * Utility Functions - Implementation
 * ================================================================ */

double lorentz_gamma(double beta)
{
    if (beta <= 0.0) return 1.0;
    if (beta >= 1.0) return INFINITY;
    return 1.0 / sqrt(1.0 - square(beta));
}

double lorentz_velocity(double px, double py, double pz, double mass)
{
    double p_mag = sqrt(square(px) + square(py) + square(pz));
    double energy = sqrt(square(px) + square(py) + square(pz) + square(mass));

    if (energy <= 0.0) return 0.0;
    return p_mag / energy;
}

double lorentz_energy(double px, double py, double pz, double mass)
{
    return sqrt(square(px) + square(py) + square(pz) + square(mass));
}

double lorentz_invariant_mass(const Particle *p1, const Particle *p2,
                               double mass)
{
    double e1, e2, px_tot, py_tot, pz_tot, s;

    if (p1 == NULL || p2 == NULL) return 0.0;

    /* Total energy */
    e1 = sqrt(square(p1->px) + square(p1->py) + square(p1->pz) + square(mass));
    e2 = sqrt(square(p2->px) + square(p2->py) + square(p2->pz) + square(mass));

    /* Total momentum */
    px_tot = p1->px + p2->px;
    py_tot = p1->py + p2->py;
    pz_tot = p1->pz + p2->pz;

    /* Invariant mass: s = (E1+E2)^2 - |p_total|^2 */
    s = square(e1 + e2) - (square(px_tot) + square(py_tot) + square(pz_tot));

    if (s < 0.0) return 0.0;  /* Avoid sqrt of negative due to rounding */
    return sqrt(s);
}

/* ================================================================
 * Lorentz Transformation Parameters
 * ================================================================ */

int lorentz_params_from_pair(const Particle *p1, const Particle *p2,
                              double mass, LorentzParams *params_out)
{
    double e1, e2, e_total;
    double px_tot, py_tot, pz_tot;
    double beta_x, beta_y, beta_z, beta_sq;
    double gamma;
    double gm1ob2, s, cm_energy;

    if (p1 == NULL || p2 == NULL || params_out == NULL) return -1;

    /* Calculate total energy */
    e1 = sqrt(square(p1->px) + square(p1->py) + square(p1->pz) + square(mass));
    e2 = sqrt(square(p2->px) + square(p2->py) + square(p2->pz) + square(mass));
    e_total = e1 + e2;

    if (e_total <= 0.0) return -1;

    /* Total momentum */
    px_tot = p1->px + p2->px;
    py_tot = p1->py + p2->py;
    pz_tot = p1->pz + p2->pz;

    /* Center-of-mass velocity: beta = p_total / E_total */
    beta_x = px_tot / e_total;
    beta_y = py_tot / e_total;
    beta_z = pz_tot / e_total;

    /* Magnitude of velocity */
    beta_sq = square(beta_x) + square(beta_y) + square(beta_z);

    /* Sanity check */
    if (beta_sq >= 1.0) return -1;

    /* Gamma factor */
    gamma = (beta_sq > 0.0) ? (1.0 / sqrt(1.0 - beta_sq)) : 1.0;

    /* Pre-compute (gamma-1)/beta^2 for use in boost formula */
    gm1ob2 = 0.0;
    if (beta_sq > 1.0e-16) {
        gm1ob2 = (gamma - 1.0) / beta_sq;
    }

    /* Invariant mass (center of mass energy) */
    s = square(e_total) - (square(px_tot) + square(py_tot) + square(pz_tot));
    if (s < 0.0) s = 0.0;
    cm_energy = sqrt(s);

    /* Fill output structure */
    params_out->cm_energy = cm_energy;
    params_out->beta_x = beta_x;
    params_out->beta_y = beta_y;
    params_out->beta_z = beta_z;
    params_out->beta = sqrt(beta_sq);
    params_out->gamma = gamma;
    params_out->gamma_m1_over_beta2 = gm1ob2;

    return 0;
}

/* ================================================================
 * Lorentz Boost to/from CM Frame
 * ================================================================
 *
 * Standard relativistic boost formula:
 *   When boosting with velocity -beta (lab to CM frame):
 *
 *   p'_mag = (gamma-1)/beta^2 * (p·beta) - gamma*E
 *   p' = p + p'_mag * beta
 *
 * This formula handles all boost directions and magnitudes correctly.
 * The factor (gamma-1)/beta^2 is numerically stable even for small beta.
 */

void lorentz_boost_to_cm(Particle *particle, const LorentzParams *params)
{
    double px, py, pz, mass;
    double energy;
    double p_dot_beta;
    double gm1ob2_p1b_mgE;
    double beta_x, beta_y, beta_z;

    if (particle == NULL || params == NULL) return;

    /* Current momentum and mass */
    px = particle->px;
    py = particle->py;
    pz = particle->pz;
    mass = NUCLEON_MASS;  /* From types.h */

    /* Total energy */
    energy = sqrt(square(px) + square(py) + square(pz) + square(mass));

    /* Components of boost velocity */
    beta_x = params->beta_x;
    beta_y = params->beta_y;
    beta_z = params->beta_z;

    /* p·beta = p_x*beta_x + p_y*beta_y + p_z*beta_z */
    p_dot_beta = px * beta_x + py * beta_y + pz * beta_z;

    /* Boost formula: p'_mag = (gamma-1)/beta^2 * (p·beta) - gamma*E */
    gm1ob2_p1b_mgE = params->gamma_m1_over_beta2 * p_dot_beta
                   - params->gamma * energy;

    /* Apply boost to momentum: p' = p + p'_mag * beta */
    particle->px = px + gm1ob2_p1b_mgE * beta_x;
    particle->py = py + gm1ob2_p1b_mgE * beta_y;
    particle->pz = pz + gm1ob2_p1b_mgE * beta_z;
}

void lorentz_boost_from_cm(Particle *particle, const LorentzParams *params)
{
    double px, py, pz, mass;
    double energy;
    double p_dot_beta;
    double gm1ob2_p1b_plus_mgE;
    double beta_x, beta_y, beta_z;

    if (particle == NULL || params == NULL) return;

    /* Current momentum and mass */
    px = particle->px;
    py = particle->py;
    pz = particle->pz;
    mass = NUCLEON_MASS;  /* From types.h */

    /* Total energy */
    energy = sqrt(square(px) + square(py) + square(pz) + square(mass));

    /* Components of boost velocity (inverse transformation uses +beta) */
    beta_x = params->beta_x;
    beta_y = params->beta_y;
    beta_z = params->beta_z;

    /* p·beta */
    p_dot_beta = px * beta_x + py * beta_y + pz * beta_z;

    /* Inverse boost formula: p'_mag = (gamma-1)/beta^2 * (p·beta) + gamma*E */
    gm1ob2_p1b_plus_mgE = params->gamma_m1_over_beta2 * p_dot_beta
                        + params->gamma * energy;

    /* Apply inverse boost to momentum: p' = p + p'_mag * beta */
    particle->px = px + gm1ob2_p1b_plus_mgE * beta_x;
    particle->py = py + gm1ob2_p1b_plus_mgE * beta_y;
    particle->pz = pz + gm1ob2_p1b_plus_mgE * beta_z;
}
