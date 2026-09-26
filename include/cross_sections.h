/*
 * LHBUU Nucleon-Nucleon Cross Sections
 *
 * Implementation: Cugnon parametrization (NIM B111(96)215)
 * Reference: J. Cugnon et al., Nucl. Instr. Meth. B 111 (1996) 215
 *
 * Physics:
 * - Energy-dependent σ_nn(√s)
 * - Isospin-dependent (pp, nn, pn channels)
 * - In-medium density reduction
 *
 * All cross sections in millibarns (mb)
 * Energies in GeV
 */

#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

#include "constants.h"
#include "macros.h"

/* Cugnon σ_nn(√s) parametrization */
double cugnon_sigma_nn(double sqrt_s, int isospin_symmetric);

/* Wrapper with density reduction */
double cugnon_sigma_nn_medium(double sqrt_s,
                               int isospin_symmetric,
                               double local_density,
                               double rho_reduction_factor);

/* Lab momentum from √s (MeV/c) */
double lab_momentum_from_sqrt_s(double sqrt_s, double nucleon_mass);

/* Isospin-averaged cross section */
double cugnon_sigma_nn_averaged(double sqrt_s);

#endif /* CROSS_SECTIONS_H */
