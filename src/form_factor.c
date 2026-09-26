/*
 * M3 Piecewise Quadratic Form Factor
 * From: J.J Monaghan, Computer Physics Reports 3 (1985) p.71
 *
 * Separable form factor: f(x,y,z) = f_x(x) * f_y(y) * f_z(z)
 * Compact support: 1.5 cells
 * C^0 continuous (continuous, but derivative has jump)
 *
 * This is a direct port from lib/form-factor.c
 */

#include "form_factor.h"
#include "constants.h"
#include "macros.h"
#include <math.h>

/* Compile-time constants - static precomputation */
static const double three_quart_dx_sq = 0.75 * DX * DX;
static const double half_dx = 0.5 * DX;
static const double three_halfs_dx = 1.5 * DX;

/*
 * 1D form factor evaluation
 *
 * Returns:
 *   > 0: form factor value
 *   -1: out of range (particle too far from grid point)
 *
 * Outputs:
 *   *sign: sign of (grid_center - particle_coord)
 *   *extras: intermediate value for gradient calculation
 */
double form_factor_1d(double particle_coord, int grid_index, double grid_center,
                      int *sign, double *extras)
{
    /* Grid cell center */
    double cell_center = (grid_index + 0.5) * DX - grid_center - particle_coord;
    double abs_dist = fabs(cell_center);

    /* Check if outside compact support (1.5 * dx) */
    double arg2 = three_halfs_dx - abs_dist;
    *extras = 0.0;
    *sign = 0;

    if (arg2 < 0.0) {
        return -1.0;  /* Out of range */
    }

    *sign = sgn(cell_center);

    /* Check if in inner region (within 0.5 * dx) */
    double arg1 = half_dx - abs_dist;

    if (arg1 < 0.0) {
        /* Outer piece: (0.5 * dx < |r| < 1.5 * dx) */
        *extras = arg2;
        return 0.5 * square(arg2);
    }

    /* Inner piece: (|r| < 0.5 * dx) */
    *extras = 2.0 * abs_dist;
    return three_quart_dx_sq - square(abs_dist);
}

/*
 * 3D form factor: product of three 1D evaluations
 *
 * Returns:
 *   > 0: form factor value f(x,y,z) = f_x * f_y * f_z
 *   -1: out of range in any dimension
 */
double form_factor_3d(double x, double y, double z,
                      int Lx, int Ly, int Lz,
                      double grid_center_x, double grid_center_y, double grid_center_z,
                      int *sign_out, double *extras_out)
{
    int sign[4];
    double extras[4];
    double g[4];  /* g[1]=x, g[2]=y, g[3]=z */

    /* Evaluate each dimension */
    g[1] = form_factor_1d(x, Lx, grid_center_x, &sign[1], &extras[1]);
    if (g[1] < 0.0) return -1.0;

    g[2] = form_factor_1d(y, Ly, grid_center_y, &sign[2], &extras[2]);
    if (g[2] < 0.0) return -1.0;

    g[3] = form_factor_1d(z, Lz, grid_center_z, &sign[3], &extras[3]);
    if (g[3] < 0.0) return -1.0;

    /* Store for gradient calculation if requested */
    if (sign_out) {
        sign_out[1] = sign[1];
        sign_out[2] = sign[2];
        sign_out[3] = sign[3];
    }

    if (extras_out) {
        extras_out[1] = extras[1];
        extras_out[2] = extras[2];
        extras_out[3] = extras[3];
    }

    /* Return product */
    return g[1] * g[2] * g[3];
}

/*
 * Gradient of 3D form factor
 *
 * Given the 1D form factors g[] and their extras[], compute:
 *   gradient[1] = df/dx = sign[1] * extras[1] * g[2] * g[3]
 *   gradient[2] = df/dy = sign[2] * extras[2] * g[3] * g[1]
 *   gradient[3] = df/dz = sign[3] * extras[3] * g[1] * g[2]
 */
void form_factor_gradient_3d(double *g, int *sign, double *extras, double *gradient)
{
    /* Gradient components using chain rule */
    gradient[1] = sgn(sign[1]) * extras[1] * g[2] * g[3];
    gradient[2] = sgn(sign[2]) * extras[2] * g[3] * g[1];
    gradient[3] = sgn(sign[3]) * extras[3] * g[1] * g[2];
}

/*
 * Normalization factor for real-space form factor
 *
 * Formula: 1 / (N_ensembles * dx^9)
 *
 * The dx^9 comes from:
 *   - dx^3 from volume element in real space
 *   - dx^3 from form factor integration
 *   - dx^3 from another form factor (for potential calculation)
 */
double form_factor_norm_r(int n_ensembles)
{
    double dx3 = cube(DX);
    double dx9 = cube(dx3);
    return 1.0 / (n_ensembles * dx9);
}
