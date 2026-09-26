/*
 * Half-Impulse Force Calculation
 *
 * Direct port from legacy/lib/half_impulse.c
 *
 * Physics: F = -∇_r U where U is total mean-field potential
 * Implementation: F = Σ_L [∇ρ(L) × U(L)] over grid cells L
 *
 * Key features from original:
 * - Self-energy correction to prevent spurious self-forces
 * - Multiple potential components (nuclear, Coulomb, isospin, momentum-dep, surface)
 * - Boundary checks for surface term
 * - Static coefficient caching for performance
 */

#include "half_impulse.h"
#include "form_factor.h"
#include "macros.h"
#include <math.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * U_coulomb: Coulomb potential at grid point
 *
 * Simple wrapper around pre-computed Coulomb grid
 * Factor of q accounts for test charge
 */
static inline double U_coulomb(const double ***coulomb_grid,
                                int l, int m, int n,
                                int q_particle)
{
    if (q_particle == 0) return 0.0;  /* Neutrons don't feel Coulomb */
    return coulomb_grid[l][m][n] * (double)q_particle;
}

/*
 * U_isospin: Isospin asymmetry potential at grid point
 *
 * From original lib/U_isospin.c
 * U_iso = D × (ρ_n - ρ_p) / ρ_0 × τ_3
 * where τ_3 = +1 for protons, -1 for neutrons
 *
 * For now: stub (requires charge density grid)
 */
static inline double U_isospin(int l, int m, int n,
                                int q_particle,
                                double D_pot)
{
    (void)l; (void)m; (void)n; (void)q_particle;
    if (D_pot == 0.0) return 0.0;
    /* TODO: Implement when charge density grid available */
    return 0.0;
}

/*
 * U_alpha_p_md: Momentum-dependent potential at grid point
 *
 * From original lib/U_alpha_p_md.c
 * Requires full momentum-space calculation (complex)
 *
 * For now: stub (requires momentum-dependent infrastructure)
 */
static inline double U_alpha_p_md(int l, int m, int n, int k,
                                   double C_pot)
{
    (void)l; (void)m; (void)n; (void)k;
    if (C_pot == 0.0) return 0.0;
    /* TODO: Implement when momentum-dependent infrastructure available */
    return 0.0;
}

/*
 * U_alpha_sur: Surface Yukawa potential at grid point
 *
 * From original lib/U_alpha_sur.c
 * Yukawa-type surface term: U_sur = A_sur × f(r) × exp(-r/λ) / r
 *
 * For now: stub (requires surface density calculation)
 */
static inline double U_alpha_sur(int l, int m, int n, int k,
                                  double A_sur)
{
    (void)l; (void)m; (void)n; (void)k;
    if (A_sur == 0.0) return 0.0;
    /* TODO: Implement when surface infrastructure available */
    return 0.0;
}

/*
 * Calculate self-energy correction
 *
 * From original lib/half_impulse.c lines 114-138
 *
 * Removes particle's own contribution to prevent self-force:
 *   correction = A × R_α(i) + σB × (ρ/ρ_0)^(σ-1) × R_α(i)
 * where R_α(i) = f_x(x_i) × f_y(y_i) × f_z(z_i) / ρ_0
 */
static double self_energy_correction(int k,
                                      int l, int m, int n,
                                      const double *x, const double *y, const double *z,
                                      double grid_center_x, double grid_center_y, double grid_center_z,
                                      const double ***rho_grid,
                                      double A_pot, double B_pot, double SIG_pot, double RHO0,
                                      double dx)
{
    (void)dx;  /* DX is compile-time constant, parameter unused */
    double correction = 0.0;
    int sign;
    double extras;

    /* Evaluate form factor at particle k's position for cell (l,m,n) */
    double R_x = form_factor_1d(x[k], l, grid_center_x, &sign, &extras);
    if (R_x <= 0) return 0.0;

    double R_y = form_factor_1d(y[k], m, grid_center_y, &sign, &extras);
    if (R_y <= 0) return 0.0;

    double R_z = form_factor_1d(z[k], n, grid_center_z, &sign, &extras);
    if (R_z <= 0) return 0.0;

    /* Form factor normalization */
    double ff_norm = form_factor_norm_r(N_ENSEMBLES);
    double R_alpha_i = R_x * R_y * R_z * ff_norm / RHO0;

    /* Skyrme self-energy: A×ρ term */
    correction = A_pot * R_alpha_i;

    /* Skyrme self-energy: B×ρ^σ term (derivative) */
    double rho = rho_grid[l][m][n];
    if (rho > 0.0) {
        double rho_ratio = rho / RHO0;
        correction += SIG_pot * B_pot * pow(rho_ratio, SIG_pot - 1.0) * R_alpha_i;
    }

    return correction;
}

/*
 * Half-impulse calculation for single particle
 *
 * Direct port from legacy/lib/half_impulse.c lines 20-163
 */
void half_impulse(int k,
                  int Lx, int Ly, int Lz,
                  double *a_part,
                  int Coul_only,
                  const DrhoDr *drho_dr,
                  const double ***U_nuc_grid,
                  const double ***coulomb_grid,
                  const double *x, const double *y, const double *z,
                  const int *q,
                  double grid_center_x, double grid_center_y, double grid_center_z,
                  int nx, int ny, int nz,
                  double dt,
                  int n_ensembles,
                  int enable_coulomb,
                  double D_pot,
                  double C_pot,
                  double A_sur,
                  const double ***rho_grid,
                  double A_pot, double B_pot, double SIG_pot, double RHO0,
                  int no_self_energy)
{
    const int ff_range = FF_RANGE;

    /* Static normalization coefficient (lines 84-89)
     * coeff = -0.5 × dt × N_ens × dx³ × form_factor_norm */
    static double half_impulse_norm = 0.0;
    if (half_impulse_norm == 0.0) {
        double coeff = -0.5 * dt * (double)n_ensembles * DX * DX * DX;
        half_impulse_norm = coeff * form_factor_norm_r(n_ensembles);
    }

    /* Search range around particle k (lines 90-95) */
    int lx_min = Lx - ff_range;
    int lx_max = Lx + ff_range;
    int ly_min = Ly - ff_range;
    int ly_max = Ly + ff_range;
    int lz_min = Lz - ff_range;
    int lz_max = Lz + ff_range;

    /* Boundary checks for surface term (lines 77-83)
     * Surface requires 2×ff_range margin */
    int sx_min = ff_range + 1;
    int sy_min = ff_range + 1;
    int sz_min = ff_range + 1;
    int sx_max = nx - ff_range;
    int sy_max = ny - ff_range;
    int sz_max = nz - ff_range;

    /* Loop over grid cells within ff_range of particle k (lines 99-159) */
    for (int l = lx_min; l <= lx_max; l++) {
        for (int m = ly_min; m <= ly_max; m++) {
            for (int n = lz_min; n <= lz_max; n++) {

                /* Get density gradient from drho_dr cache (lines 103-108)
                 * gradient[0] = ∂f/∂x, gradient[1] = ∂f/∂y, gradient[2] = ∂f/∂z
                 * Skip if gradient is zero (particle doesn't contribute to this cell) */
                double gradient[3];
                gradient[0] = drho_dr_get(drho_dr, k, l - lx_min, m - ly_min, n - lz_min, 0);
                if (gradient[0] == 0.0) continue;

                gradient[1] = drho_dr_get(drho_dr, k, l - lx_min, m - ly_min, n - lz_min, 1);
                if (gradient[1] == 0.0) continue;

                gradient[2] = drho_dr_get(drho_dr, k, l - lx_min, m - ly_min, n - lz_min, 2);
                if (gradient[2] == 0.0) continue;

                /* Calculate total single-particle potential U(l,m,n) (lines 109-155) */
                double U_single = 0.0;

                if (!Coul_only) {
                    /* Nuclear mean-field (line 112) */
                    U_single = U_nuc_grid[l][m][n];

                    /* Self-energy correction (lines 114-138) */
                    if (no_self_energy == 1) {
                        double correction = self_energy_correction(k, l, m, n,
                                                                   x, y, z,
                                                                   grid_center_x, grid_center_y, grid_center_z,
                                                                   rho_grid,
                                                                   A_pot, B_pot, SIG_pot, RHO0,
                                                                   DX);
                        U_single -= correction;
                    }

                    /* Coulomb potential (line 140) */
                    if (enable_coulomb && q[k]) {
                        U_single += U_coulomb(coulomb_grid, l, m, n, q[k]);
                    }

                    /* Isospin asymmetry potential (line 142) */
                    if (D_pot != 0.0) {
                        U_single += U_isospin(l, m, n, q[k], D_pot);
                    }

                    /* Momentum-dependent potential (line 144) */
                    if (C_pot != 0.0) {
                        U_single += U_alpha_p_md(l, m, n, k, C_pot);
                    }

                    /* Surface Yukawa potential (lines 146-150) */
                    if (A_sur != 0.0) {
                        if (sx_min <= l && l <= sx_max &&
                            sy_min <= m && m <= sy_max &&
                            sz_min <= n && n <= sz_max) {
                            U_single += U_alpha_sur(l, m, n, k, A_sur);
                        }
                    }
                } else {
                    /* Coulomb-only mode (line 154) */
                    if (enable_coulomb && q[k]) {
                        U_single = U_coulomb(coulomb_grid, l, m, n, q[k]);
                    }
                }

                /* Accumulate force: F += ∇ρ × U (lines 156-158) */
                a_part[0] += gradient[0] * U_single;
                a_part[1] += gradient[1] * U_single;
                a_part[2] += gradient[2] * U_single;
            }
        }
    }

    /* Apply normalization (lines 160-162) */
    a_part[0] *= half_impulse_norm;
    a_part[1] *= half_impulse_norm;
    a_part[2] *= half_impulse_norm;
}

/*
 * Batch calculation: all particles with OpenMP parallelization
 *
 * OpenMP strategy:
 * - Parallel over particles (outer loop)
 * - Dynamic scheduling for load balance
 * - Chunk size 64 balances overhead vs. granularity
 */
void half_impulse_all(int n_particles,
                      const int *Lx_arr, const int *Ly_arr, const int *Lz_arr,
                      double **a_part_arr,
                      int Coul_only,
                      const DrhoDr *drho_dr,
                      const double ***U_nuc_grid,
                      const double ***coulomb_grid,
                      const double *x, const double *y, const double *z,
                      const int *q,
                      double grid_center_x, double grid_center_y, double grid_center_z,
                      int nx, int ny, int nz,
                      double dt,
                      int n_ensembles,
                      int enable_coulomb,
                      double D_pot,
                      double C_pot,
                      double A_sur,
                      const double ***rho_grid,
                      double A_pot, double B_pot, double SIG_pot, double RHO0,
                      int no_self_energy)
{
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
#endif
    for (int k = 1; k <= n_particles; k++) {
        /* Pre-zero acceleration (required by half_impulse) */
        a_part_arr[k][0] = 0.0;
        a_part_arr[k][1] = 0.0;
        a_part_arr[k][2] = 0.0;

        /* Calculate force */
        half_impulse(k,
                     Lx_arr[k], Ly_arr[k], Lz_arr[k],
                     a_part_arr[k],
                     Coul_only,
                     drho_dr,
                     U_nuc_grid,
                     coulomb_grid,
                     x, y, z,
                     q,
                     grid_center_x, grid_center_y, grid_center_z,
                     nx, ny, nz,
                     dt,
                     n_ensembles,
                     enable_coulomb,
                     D_pot,
                     C_pot,
                     A_sur,
                     rho_grid,
                     A_pot, B_pot, SIG_pot, RHO0,
                     no_self_energy);
    }
}
