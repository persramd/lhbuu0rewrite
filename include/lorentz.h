/*
 * LHBUU Lorentz Transformation Utilities
 *
 * Provides boost transformations to/from center-of-mass frame
 * for nucleon-nucleon collision physics.
 *
 * Based on original lorentz_boostnntoNN.c and collisions.c
 * from /Users/declan/Documents/dev/osx/lhbuu/lib/
 *
 * All calculations use double precision for numerical stability.
 */

#ifndef LHBUU_LORENTZ_H
#define LHBUU_LORENTZ_H

#include "types.h"

/* ================================================================
 * Lorentz Parameters Structure
 * ================================================================
 *
 * Encapsulates the transformation parameters needed to boost
 * between lab frame and center-of-mass frame.
 *
 * Fields:
 *   cm_energy     : Invariant mass of the system (sqrt(s))
 *   beta_x, y, z  : Velocity components of CM frame (units of c)
 *   beta           : Magnitude of velocity |beta|
 *   gamma          : Lorentz factor 1/sqrt(1-beta^2)
 */
typedef struct {
    double cm_energy;
    double beta_x;
    double beta_y;
    double beta_z;
    double beta;
    double gamma;
    double gamma_m1_over_beta2;  /* (gamma-1)/beta^2 for boost formula */
} LorentzParams;

/* ================================================================
 * Core Lorentz Transformations
 * ================================================================ */

/*
 * Calculate center-of-mass parameters from two particles.
 *
 * Computes the CM frame velocity and Lorentz factor from the
 * momenta of two colliding nucleons.
 *
 * Args:
 *   p1, p2     : Particles with px, py, pz, mass
 *   mass       : Nucleon mass (MeV/c²)
 *   params_out : Output structure with CM parameters
 *
 * Returns: 0 on success, -1 on error (beta >= 1)
 */
int lorentz_params_from_pair(const Particle *p1, const Particle *p2,
                              double mass, LorentzParams *params_out);

/*
 * Boost a particle's 4-momentum to the center-of-mass frame.
 *
 * Applies Lorentz transformation to a particle in the lab frame,
 * moving it to the rest frame of the two-particle system.
 *
 * Formula (standard relativistic boost):
 *   p'_parallel = (gamma-1)/beta^2 * (p·beta) - gamma*E
 *   p' = p + p'_parallel * beta
 *
 * Args:
 *   particle   : Particle to transform (modified in-place)
 *   params     : Lorentz parameters of the boost
 */
void lorentz_boost_to_cm(Particle *particle, const LorentzParams *params);

/*
 * Boost a particle's 4-momentum from CM frame back to lab frame.
 *
 * Inverse transformation: applies boost with -beta direction.
 * Used after scattering calculation to return particles to lab frame.
 *
 * Args:
 *   particle   : Particle to transform (modified in-place)
 *   params     : Lorentz parameters of the original boost
 */
void lorentz_boost_from_cm(Particle *particle, const LorentzParams *params);

/* ================================================================
 * Utility Functions
 * ================================================================ */

/*
 * Calculate the Lorentz factor gamma.
 *
 * gamma = 1 / sqrt(1 - beta^2)
 *
 * Args:
 *   beta : Velocity magnitude (0 < beta < 1)
 *
 * Returns: Lorentz factor (always >= 1.0)
 */
double lorentz_gamma(double beta);

/*
 * Calculate velocity magnitude from momentum and mass.
 *
 * beta = |p| / E = |p| / sqrt(|p|^2 + m^2)
 *
 * Args:
 *   px, py, pz : Momentum components (MeV/c)
 *   mass       : Rest mass (MeV/c²)
 *
 * Returns: Velocity magnitude in units of c (0 <= beta < 1)
 */
double lorentz_velocity(double px, double py, double pz, double mass);

/*
 * Calculate total energy from momentum and mass.
 *
 * E = sqrt(|p|^2 + m^2)
 *
 * Args:
 *   px, py, pz : Momentum components (MeV/c)
 *   mass       : Rest mass (MeV/c²)
 *
 * Returns: Total energy (MeV)
 */
double lorentz_energy(double px, double py, double pz, double mass);

/*
 * Calculate invariant mass (sqrt(s)) of a two-particle system.
 *
 * The invariant mass is frame-independent and equals the CM energy.
 *
 * s = (E1 + E2)^2 - |p1 + p2|^2
 * sqrt(s) = CM energy
 *
 * Args:
 *   p1, p2  : Particles with px, py, pz, mass
 *   mass    : Nucleon mass (MeV/c²)
 *
 * Returns: Invariant mass (MeV)
 */
double lorentz_invariant_mass(const Particle *p1, const Particle *p2,
                               double mass);

#endif /* LHBUU_LORENTZ_H */
