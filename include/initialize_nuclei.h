/*
 * LHBUU Woods-Saxon Nuclear Initialization
 *
 * Implements Woods-Saxon density profile for realistic nuclear initialization
 * Ported from original LHBUU lib/r_p_initialization.c
 *
 * Key features:
 * - Woods-Saxon spatial density: rho(r) = rho0 / (1 + exp((r-R)/a))
 * - Local Thomas-Fermi momentum sampling with Lenk-Pandharipande corrections
 * - Ensemble generation for statistical averaging
 * - Double precision for numerical stability
 */

#ifndef INITIALIZE_NUCLEI_H
#define INITIALIZE_NUCLEI_H

#include "types.h"
#include "constants.h"
#include "macros.h"

/* Woods-Saxon parameters */
typedef struct {
    double R;              /* Nuclear radius (fm) */
    double a;              /* Surface diffuseness (fm) */
    double rho0;           /* Central density (fm^-3) */
    double skin_thickness; /* Skin thickness parameter (fm) */
} WoodsSaxonParams;

/* Nucleus initialization configuration */
typedef struct {
    int A;                 /* Mass number */
    int Z;                 /* Atomic number */
    int N_ensembles;       /* Number of ensembles */
    WoodsSaxonParams ws;   /* Woods-Saxon parameters */
    double min_separation; /* Minimum particle separation (fm) */
    int is_projectile;     /* 1=projectile, 0=target */
} NucleusInitConfig;

/* Function prototypes */

/*
 * Initialize Woods-Saxon parameters for a nucleus
 * Uses standard empirical formulae:
 *   R = 1.12 * A^(1/3) fm
 *   a = 0.54 fm (surface diffuseness)
 *   rho0 = 0.16 fm^-3 (nuclear saturation density)
 */
void init_woods_saxon_params(WoodsSaxonParams *ws, int A);

/*
 * Sample radius from Woods-Saxon density profile
 * Uses rejection sampling with normalized profile
 * Returns: radial distance in fm
 */
double sample_woods_saxon_radius(const WoodsSaxonParams *ws);

/*
 * Calculate Woods-Saxon density at radius r
 * rho(r) = rho0 / (1 + exp((r-R)/a))
 */
double woods_saxon_density(double r, const WoodsSaxonParams *ws);

/*
 * Initialize nuclear positions using Woods-Saxon profile
 * Ensures minimum separation between particles
 * Returns: number of particles initialized
 */
int initialize_nuclear_positions(
    Particle *particles,
    const NucleusInitConfig *config,
    int start_idx
);

/*
 * Initialize nuclear momenta using Local Thomas-Fermi approximation
 * Implements Lenk-Pandharipande corrections (PRC39, 2242, 1989)
 * Requires: positions already initialized, density grid computed
 * Returns: 0 on success, -1 on error
 */
int initialize_nuclear_momenta(
    Particle *particles,
    int num_particles,
    const Grid *grid,
    double D_pot  /* Isospin potential strength */
);

/*
 * Calculate local Fermi momentum with LTF corrections
 * Implements quadratic spline form factor (A_ff = 17/288)
 * Returns: local Fermi momentum in MeV
 */
double calculate_local_fermi_momentum(
    double rho_local,
    double rho_nabla2,   /* Laplacian term */
    double rho_grad2,    /* Gradient squared term */
    int is_charged,      /* 1=proton, 0=neutron */
    double D_pot
);

/*
 * Zero center-of-mass momentum for nucleus
 * Ensures overall momentum conservation
 */
void zero_cm_momentum(
    Particle *particles,
    int start_idx,
    int count,
    int N_ensembles
);

/*
 * Generate full nucleus (positions + momenta)
 * High-level initialization routine
 * Returns: 0 on success, -1 on error
 */
int generate_nucleus(
    Particle *particles,
    const NucleusInitConfig *config,
    const Grid *grid,
    int start_idx,
    double D_pot
);

#endif /* INITIALIZE_NUCLEI_H */
