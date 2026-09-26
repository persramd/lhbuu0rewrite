/*
 * LHBUU Compile-Time Constants
 *
 * Philosophy: Recompile for different parameters
 * - Enables aggressive compiler optimization
 * - Static precomputation at compile time
 * - Architecture-specific tuning
 *
 * To change parameters: edit this file and recompile
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>

/* =================================================================
 * Grid Parameters - MUST RECOMPILE FOR DIFFERENT SYSTEMS
 * ================================================================= */

/* Grid spacing (fm) */
#define DX 1.0

/* Grid dimensions (number of cells) */
#define NX 64
#define NY 64
#define NZ 64

/* Total grid cells */
#define N_CELLS (NX * NY * NZ)

/* Form factor range (±2 cells) */
#define FF_RANGE 2

/* Particles per grid cell (for allocation) */
#define PARTICLES_PER_GRID 16

/* =================================================================
 * Physical Constants
 * ================================================================= */

/* Note: NUCLEON_MASS and HBARC defined in types.h (GeV units) */

/* Nucleon mass (MeV/c²) - for legacy code compatibility */
#define NUCLEON_MASS_MEV 938.0

/* Pion mass (MeV/c²) */
#define PION_MASS 138.0

/* Fine structure constant (dimensionless) */
#define ALPHA_EM (1.0/137.0)

/* Hbar*c (MeV·fm) - for legacy code compatibility */
#define HBARC_MEV 197.33

/* Nuclear saturation density (fm⁻³) */
#define RHO_0 0.16

/* Fermi momentum at saturation (MeV/c) */
#define P_F0_MEV 268.0

/* =================================================================
 * Potential Parameters
 * ================================================================= */

/* Skyrme parameters (MeV) */
#define SKYRME_A -356.0
#define SKYRME_B 303.0
#define SKYRME_SIGMA 7.0/6.0
#define SKYRME_TAU 2.0

/* Momentum-dependent potential parameters */
#define C_MDYI 2.0        /* Strength (dimensionless) */
#define LAMBDA_MDYI 1.5   /* Range parameter */

/* Isospin asymmetry energy (MeV) */
#define C_SYM 32.0

/* =================================================================
 * Collision Parameters
 * ================================================================= */

/* Cugnon cross section parameters */
#define SIGMA_0 40.0      /* Free NN cross section (mb) */

/* In-medium reduction factor */
#define ALPHA_MEDIUM 0.2

/* Pauli blocking parameters */
#define PAULI_RADIUS 2.0  /* Spatial radius (fm) */
#define PAULI_P_RADIUS (200.0/HBARC)  /* Momentum radius (fm⁻¹) */

/* =================================================================
 * Ensemble Parameters
 * ================================================================= */

/* Number of ensembles */
#define N_ENSEMBLES 100

/* Maximum particles per ensemble */
#define MAX_PARTICLES_PER_ENSEMBLE 200

/* Total maximum particles */
#define MAX_PARTICLES (N_ENSEMBLES * MAX_PARTICLES_PER_ENSEMBLE)

/* =================================================================
 * Integration Parameters
 * ================================================================= */

/* Time step (fm/c) */
#define DT 0.2

/* Maximum number of time steps */
#define MAX_STEPS 1000

/* Diagnostic output interval */
#define DIAG_INTERVAL 10

/* Grid sort interval (0 = never) */
#define SORT_INTERVAL 50

/* =================================================================
 * Derived Constants (DO NOT EDIT - computed at compile time)
 * ================================================================= */

/* Grid volume */
#define GRID_VOLUME (DX * DX * DX)

/* Inverse grid spacing */
#define INV_DX (1.0/DX)

/* Lambda squared for momentum-dependent potential */
#define LAMBDA_SQ (LAMBDA_MDYI * LAMBDA_MDYI * P_F0_MEV * P_F0_MEV)

/* Form factor normalization */
#define FF_NORM (1.0/(N_ENSEMBLES * GRID_VOLUME * GRID_VOLUME * GRID_VOLUME))

#endif /* CONSTANTS_H */
