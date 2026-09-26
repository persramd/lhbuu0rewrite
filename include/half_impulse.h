/*
 * Half-Impulse Force Calculation
 *
 * Physics:
 * --------
 * Computes mean-field force F = -∇_r U on each particle from the gradient
 * of the total mean-field potential U(r). The "half-impulse" name comes from
 * the leapfrog integrator where forces are evaluated at half-timestep intervals.
 *
 * Force components:
 *   F_i = -∇_r [U_nuclear + U_coulomb + U_isospin + U_momentum_dep + U_surface]
 *
 * Algorithm (from legacy lib/half_impulse.c):
 * --------------------------------------------
 * 1. Loop over grid cells L within ff_range of particle k
 * 2. For each cell, get pre-computed density gradient: ∇ρ_k(L) from drho_dr
 * 3. Evaluate total potential U(L) including self-energy correction
 * 4. Accumulate force: F += ∇ρ_k(L) × U(L)
 * 5. Apply normalization: F *= -0.5 × dt × N_ens × dx³ × form_factor_norm
 *
 * Self-energy correction (NO_SELF_ENERGY flag):
 *   Subtract particle's own contribution to prevent spurious self-forces
 *
 * Port status: Direct from legacy/lib/half_impulse.c (lines 20-163)
 */

#ifndef HALF_IMPULSE_H
#define HALF_IMPULSE_H

#include "constants.h"
#include "drho_dr.h"

/*
 * Calculate half-impulse force on single particle
 *
 * From original: lib/half_impulse.c lines 20-163
 *
 * Args:
 *   k: particle index (1-based)
 *   Lx, Ly, Lz: particle's grid cell location
 *   a_part: (output) acceleration [0..2] (MUST be pre-zeroed by caller)
 *   Coul_only: 1=Coulomb only, 0=all potentials
 *   drho_dr: pre-computed density gradients
 *   U_nuc_grid: nuclear mean-field potential on grid [nx][ny][nz]
 *   coulomb_grid: Coulomb potential on grid [nx][ny][nz]
 *   x, y, z: particle positions [1..n_particles]
 *   q: particle charges [1..n_particles] (0=neutron, 1=proton)
 *   grid_center_x, grid_center_y, grid_center_z: grid origins
 *   nx, ny, nz: grid dimensions
 *   dt: timestep (fm/c)
 *   n_ensembles: number of ensembles
 *   enable_coulomb: 1=include Coulomb, 0=skip
 *   D_pot: isospin asymmetry strength (MeV)
 *   C_pot: momentum-dependent strength
 *   A_sur: surface Yukawa strength (MeV)
 *   rho_grid: baryon density on grid [nx][ny][nz] (for self-energy)
 *   A_pot, B_pot, SIG_pot: Skyrme parameters (for self-energy)
 *   RHO0: saturation density (for self-energy)
 *   no_self_energy: 1=apply self-energy correction, 0=skip
 *
 * Output:
 *   a_part[0..2]: force components (CALLER MUST PRE-ZERO)
 *
 * Boundary handling:
 *   Returns zero force if particle too close to grid boundary (< ff_range)
 *
 * OpenMP parallelization:
 *   Outer loop over particles is parallelized
 *   Inner loops over cells are sequential (cache coherence)
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
                  int no_self_energy);

/*
 * Batch calculation: all particles
 *
 * Wrapper for parallel execution over all particles
 * Handles OpenMP threading and ensures proper data partitioning
 *
 * Args:
 *   n_particles: total number of particles
 *   Lx, Ly, Lz: particle grid locations [1..n_particles]
 *   a_part: (output) accelerations [particle][0..2] (MUST be pre-zeroed)
 *   ... (same parameters as half_impulse)
 *
 * OpenMP strategy:
 *   #pragma omp parallel for schedule(dynamic, 64)
 *   Dynamic scheduling for load balance (particles have varying neighbors)
 *   Chunk size 64 balances overhead vs. granularity
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
                      int no_self_energy);

#endif /* HALF_IMPULSE_H */
