/*
 * LHBUU Nucleon-Nucleon Cross Sections
 *
 * Implementation: Cugnon parametrization (NIM B111(96)215)
 *
 * Original implementation: scatter.c lines 74-120
 *
 * Physics notes:
 * - Four energy regimes with different functional forms
 * - Separate parametrizations for pp/nn (isospin symmetric) vs pn
 * - Low energy: constant σ = 150 mb
 * - Intermediate: polynomial in lab momentum
 * - High energy: 77/(p_lab + 1.5)
 *
 * Units:
 * - sqrt_s: GeV (invariant mass)
 * - p_lab: GeV/c (lab momentum)
 * - σ_nn: mb (millibarns)
 */

#include <math.h>
#include "cross_sections.h"
#include "constants.h"
#include "macros.h"

/* Threshold energies for Cugnon parametrization (GeV) */
#define SQRT_S_THRESH_1_SYM  1.8835  /* pp/nn threshold 1 */
#define SQRT_S_THRESH_1_ASYM 1.8877  /* pn threshold 1 */
#define SQRT_S_THRESH_2      2.0180  /* threshold 2 (both channels) */
#define SQRT_S_THRESH_3      2.4298  /* threshold 3 (both channels) */

/* Constant cross section at low energy (mb) */
#define SIGMA_LOW_ENERGY 150.0

/*
 * Calculate lab momentum from invariant mass √s
 *
 * For equal-mass particles:
 * p_lab = M * sqrt((s/(2M²) - 1)² - 1)
 *
 * where M is the nucleon mass.
 *
 * Derivation: scatter.c lines 81-82
 */
double lab_momentum_from_sqrt_s(double sqrt_s, double nucleon_mass)
{
    double s_over_2m2 = square(sqrt_s) / (2.0 * square(nucleon_mass));
    double term = square(s_over_2m2 - 1.0) - 1.0;

    if (term < 0.0) return 0.0;  /* Below threshold */

    return nucleon_mass * sqrt(term);
}

/*
 * Cugnon σ_nn parametrization (free space)
 *
 * Parameters:
 *   sqrt_s: Center-of-mass energy (GeV)
 *   isospin_symmetric: 1 for pp/nn, 0 for pn
 *
 * Returns: Cross section in millibarns (mb)
 *
 * Original: scatter.c lines 83-120
 */
double cugnon_sigma_nn(double sqrt_s, int isospin_symmetric)
{
    double sigma;
    double p_lab;
    const double M_NUCLEON = 0.938919;  /* GeV (from original N_MASS) */

    /* Convert √s to lab momentum (GeV/c) */
    p_lab = lab_momentum_from_sqrt_s(sqrt_s, M_NUCLEON);

    /* Isospin-symmetric channel (pp or nn) */
    if (isospin_symmetric) {
        if (sqrt_s <= SQRT_S_THRESH_1_SYM) {
            /* Low energy: constant */
            sigma = SIGMA_LOW_ENERGY;
        }
        else if (sqrt_s <= SQRT_S_THRESH_2) {
            /* Regime 2: polynomial rise */
            sigma = 23.5 + 1000.0 * pow(p_lab - 0.7, 4.0);
        }
        else if (sqrt_s <= SQRT_S_THRESH_3) {
            /* Regime 3: peak and fall */
            sigma = 1250.0 / (p_lab + 50.0) - 4.0 * square(p_lab - 1.3);
        }
        else {
            /* High energy: power law */
            sigma = 77.0 / (p_lab + 1.5);
        }
    }
    /* Isospin-asymmetric channel (pn) */
    else {
        if (sqrt_s <= SQRT_S_THRESH_1_ASYM) {
            /* Low energy: constant */
            sigma = SIGMA_LOW_ENERGY;
        }
        else if (sqrt_s <= SQRT_S_THRESH_2) {
            /* Regime 2: different polynomial */
            sigma = 33.0 + 196.0 * pow(fabs(p_lab - 0.95), 2.5);
        }
        else if (sqrt_s <= SQRT_S_THRESH_3) {
            /* Regime 3: inverse square root */
            sigma = 31.0 / sqrt(p_lab);
        }
        else {
            /* High energy: same as symmetric */
            sigma = 77.0 / (p_lab + 1.5);
        }
    }

    return sigma;
}

/*
 * Isospin-averaged cross section
 *
 * Used when isospin dependence is disabled
 *
 * σ_avg = (σ_pp + σ_pn) / 2
 *
 * Original: scatter.c lines 106-120
 */
double cugnon_sigma_nn_averaged(double sqrt_s)
{
    double sigma;
    double p_lab;
    const double M_NUCLEON = 0.938919;  /* GeV (from original N_MASS) */

    /* Convert √s to lab momentum */
    p_lab = lab_momentum_from_sqrt_s(sqrt_s, M_NUCLEON);

    /* Low energy: use one threshold (arbitrary choice in original) */
    if (sqrt_s <= SQRT_S_THRESH_1_ASYM) {
        sigma = SIGMA_LOW_ENERGY;
    }
    /* Regime 2: average both channels */
    else if (sqrt_s <= SQRT_S_THRESH_2) {
        double sigma_sym = 23.5 + 1000.0 * pow(p_lab - 0.7, 4.0);
        double sigma_asym = 33.0 + 196.0 * pow(fabs(p_lab - 0.95), 2.5);
        sigma = (sigma_sym + sigma_asym) / 2.0;
    }
    /* Regime 3: average both channels */
    else if (sqrt_s <= SQRT_S_THRESH_3) {
        double sigma_sym = 1250.0 / (p_lab + 50.0) - 4.0 * square(p_lab - 1.3);
        double sigma_asym = 31.0 / sqrt(p_lab);
        sigma = (sigma_sym + sigma_asym) / 2.0;
    }
    /* High energy: identical in both channels */
    else {
        sigma = 77.0 / (p_lab + 1.5);
    }

    return sigma;
}

/*
 * In-medium cross section with density reduction
 *
 * σ_medium = σ_free * (1 - α * ρ/ρ₀)
 *
 * Parameters:
 *   sqrt_s: Center-of-mass energy (GeV)
 *   isospin_symmetric: 1 for pp/nn, 0 for pn
 *   local_density: Local density (fm^-3)
 *   rho_reduction_factor: α parameter (typically 0.2)
 *
 * Original: scatter.c lines 134-135
 */
double cugnon_sigma_nn_medium(double sqrt_s,
                               int isospin_symmetric,
                               double local_density,
                               double rho_reduction_factor)
{
    double sigma_free;
    double sigma_medium;
    double reduction;

    /* Get free-space cross section */
    sigma_free = cugnon_sigma_nn(sqrt_s, isospin_symmetric);

    /* Apply density reduction */
    reduction = 1.0 - rho_reduction_factor * local_density;
    sigma_medium = sigma_free * reduction;

    /* Prevent negative cross sections */
    if (sigma_medium < 0.0) sigma_medium = 0.0;

    return sigma_medium;
}
