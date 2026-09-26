/*
 * M3 Piecewise Quadratic Form Factor
 *
 * Separable, compact support (1.5 cells), C^0 continuous
 * Direct port from original lib/form-factor.c
 */

#ifndef FORM_FACTOR_H
#define FORM_FACTOR_H

/*
 * 1D form factor evaluation
 *
 * Args:
 *   particle_coord: particle position in this dimension
 *   grid_index: grid cell index (L)
 *   grid_center: grid origin
 *   sign: (output) sign for gradient calculation
 *   extras: (output) intermediate value for gradient
 *
 * Returns:
 *   > 0: form factor value
 *   -1: out of range
 *
 * Note: DX is compile-time constant
 */
double form_factor_1d(double particle_coord, int grid_index, double grid_center,
                      int *sign, double *extras);

/*
 * 3D form factor: f(x,y,z) = f_x(x) * f_y(y) * f_z(z)
 *
 * Args:
 *   x, y, z: particle position
 *   Lx, Ly, Lz: grid cell indices
 *   grid_center_x, grid_center_y, grid_center_z: grid origins
 *   sign_out: (output) signs for gradient [1..3]
 *   extras_out: (output) extras for gradient [1..3]
 *
 * Returns:
 *   > 0: form factor value
 *   -1: out of range in any dimension
 *
 * Note: DX is compile-time constant
 */
double form_factor_3d(double x, double y, double z,
                      int Lx, int Ly, int Lz,
                      double grid_center_x, double grid_center_y, double grid_center_z,
                      int *sign_out, double *extras_out);

/*
 * Gradient of 3D form factor
 *
 * Args:
 *   g: 1D form factor values [1..3]
 *   sign: signs from 1D evaluation [1..3]
 *   extras: extras from 1D evaluation [1..3]
 *   gradient: (output) gradient components [1..3]
 */
void form_factor_gradient_3d(double *g, int *sign, double *extras, double *gradient);

/*
 * Normalization factor for real-space form factor
 *
 * Args:
 *   n_ensembles: number of ensembles
 *
 * Returns:
 *   1 / (N_ensembles * dx^9)
 *
 * Note: DX is compile-time constant
 */
double form_factor_norm_r(int n_ensembles);

#endif /* FORM_FACTOR_H */
