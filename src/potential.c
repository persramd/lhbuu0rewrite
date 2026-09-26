#include "potential.h"
#include "grid.h"
#include <math.h>

void potential_calculate_meanfield(Grid *grid, const PotentialParams *params) {
    size_t n_cells = grid->nx * grid->ny * grid->nz;

    for (size_t i = 0; i < n_cells; i++) {
        float rho = grid->baryon_density[i];
        /* Skyrme potential: U = A*rho + B*rho^sigma */
        grid->potential[i] = params->A * rho + params->B * powf(rho, params->sigma);
    }
}

float potential_at_point(const Grid *grid, const PotentialParams *params,
                         int ix, int iy, int iz) {
    if (ix < 0 || ix >= grid->nx || iy < 0 || iy >= grid->ny || iz < 0 || iz >= grid->nz) {
        return 0.0f;
    }
    int idx = grid_index(grid, ix, iy, iz);
    return grid->potential[idx];
}

void potential_calculate_force(const Grid *grid, const PotentialParams *params,
                                const Particle *p, float *fx, float *fy, float *fz) {
    int ix, iy, iz;
    grid_coords(grid, p->x, p->y, p->z, &ix, &iy, &iz);

    if (ix < 1 || ix >= grid->nx-1 || iy < 1 || iy >= grid->ny-1 || iz < 1 || iz >= grid->nz-1) {
        *fx = *fy = *fz = 0.0f;
        return;
    }

    /* Finite difference gradient */
    float U_xp = potential_at_point(grid, params, ix+1, iy, iz);
    float U_xm = potential_at_point(grid, params, ix-1, iy, iz);
    float U_yp = potential_at_point(grid, params, ix, iy+1, iz);
    float U_ym = potential_at_point(grid, params, ix, iy-1, iz);
    float U_zp = potential_at_point(grid, params, ix, iy, iz+1);
    float U_zm = potential_at_point(grid, params, ix, iy, iz-1);

    *fx = -(U_xp - U_xm) / (2.0f * grid->dx);
    *fy = -(U_yp - U_ym) / (2.0f * grid->dx);
    *fz = -(U_zp - U_zm) / (2.0f * grid->dx);

    /* Add Coulomb force if enabled */
    if (params->enable_coulomb && p->charge) {
        int idx = grid_index(grid, ix, iy, iz);
        int idx_xp = grid_index(grid, ix+1, iy, iz);
        int idx_xm = grid_index(grid, ix-1, iy, iz);
        int idx_yp = grid_index(grid, ix, iy+1, iz);
        int idx_ym = grid_index(grid, ix, iy-1, iz);
        int idx_zp = grid_index(grid, ix, iy, iz+1);
        int idx_zm = grid_index(grid, ix, iy, iz-1);

        float Vc_xp = grid->coulomb_potential[idx_xp];
        float Vc_xm = grid->coulomb_potential[idx_xm];
        float Vc_yp = grid->coulomb_potential[idx_yp];
        float Vc_ym = grid->coulomb_potential[idx_ym];
        float Vc_zp = grid->coulomb_potential[idx_zp];
        float Vc_zm = grid->coulomb_potential[idx_zm];

        *fx -= (Vc_xp - Vc_xm) / (2.0f * grid->dx);
        *fy -= (Vc_yp - Vc_ym) / (2.0f * grid->dx);
        *fz -= (Vc_zp - Vc_zm) / (2.0f * grid->dx);
    }
}

/* Momentum-dependent potential gradient: ∇_p V(r,p)
 * From legacy: grad_V_alpha_p_md() in lib/U_alpha_p_md.c
 *
 * This computes the gradient with respect to MOMENTUM of the
 * momentum-dependent interaction potential. This term contributes
 * to the POSITION update in the modified Verlet algorithm:
 *
 *   r(t+dt) = r(t) + v*dt + ∇_p V*dt
 *
 * where ∇_p V = ∂V/∂p for momentum-dependent potential V(r,p).
 *
 * TODO: Implement full MDYI momentum gradient calculation
 * Currently returns zero (no momentum-dependent contribution)
 */
void potential_grad_momentum_dep(const Grid *grid, const PotentialParams *params,
                                 const Particle *p, int ix, int iy, int iz,
                                 float *gx, float *gy, float *gz) {
    /* Initialize to zero */
    *gx = *gy = *gz = 0.0f;

    /* Return early if momentum-dependent potential disabled */
    if (!params->enable_momentum_dep || params->C == 0.0f) {
        return;
    }

    /* TODO: Full implementation requires:
     * 1. Loop over neighboring particles (via grid_tag or cell lists)
     * 2. For each pair (i,j): compute form factors f(r_ij)
     * 3. Calculate g2 = 1 / (λ²p_f0² + |p_i - p_j|²)
     * 4. Sum gradient: ∇_p V ∝ Σ_j f(r_ij) * g2 * (p_i - p_j)
     * 5. Apply normalization constants
     *
     * See legacy lib/U_alpha_p_md.c:344-516 for reference
     * Uses cached local_pig and g_store arrays for performance
     */
}

float potential_momentum_dependent(const Grid *grid, const PotentialParams *params,
                                   const Particle *p, int ix, int iy, int iz) {
    /* TODO: Implement MDYI potential U(r,p)
     * See legacy lib/U_alpha_p_md.c:13-161
     * Returns momentum-dependent contribution to potential energy
     */
    return 0.0f;
}
