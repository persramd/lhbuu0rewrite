/*
 * LHBUU Woods-Saxon Nuclear Initialization Implementation
 *
 * Ported from original LHBUU:
 *   - lib/r_p_initialization.c (spatial initialization)
 *   - lib/resample_rogue.c (momentum resampling)
 *
 * References:
 *   - Lenk & Pandharipande, PRC39, 2242 (1989) - LTF prescription
 *   - Cugnon et al., NPA 352, 505 (1981) - BUU formalism
 */

#include "initialize_nuclei.h"
#include "grid.h"
#include "utils.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* LTF form factor coefficients (from Lenk-Pandharipande) */
#define LTF_A2 (1.0/54.0)        /* Linear spline */
#define LTF_A3 (17.0/288.0)      /* Quadratic spline - used here */
#define LTF_A4 0.0               /* Cubic spline (not computed) */
#define LTF_A5 0.0               /* Quartic spline (not computed) */
#define A_FF LTF_A3              /* Match with form factor */

/* Minimum density threshold (1% of rho0) */
#define RHO_MIN_THRESHOLD 1e-3

/* Maximum initialization attempts */
#define MAX_INIT_ATTEMPTS 10000

/*
 * Initialize Woods-Saxon parameters
 */
void init_woods_saxon_params(WoodsSaxonParams *ws, int A) {
    /* Standard Woods-Saxon parameterization */
    ws->R = 1.12 * pow((double)A, 1.0/3.0);  /* Nuclear radius (fm) */
    ws->a = 0.54;                             /* Surface diffuseness (fm) */
    ws->rho0 = RHO_0;                         /* Saturation density */
    ws->skin_thickness = 2.3;                 /* Skin thickness (fm) */
}

/*
 * Calculate Woods-Saxon density
 * rho(r) = rho0 / (1 + exp((r-R)/a))
 */
double woods_saxon_density(double r, const WoodsSaxonParams *ws) {
    return ws->rho0 / (1.0 + exp((r - ws->R) / ws->a));
}

/*
 * Sample radius from Woods-Saxon distribution
 * Uses theta-function approximation (uniform sampling in sphere)
 * More sophisticated rejection sampling could be implemented if needed
 */
double sample_woods_saxon_radius(const WoodsSaxonParams *ws) {
    /* Simple theta-function initialization (as in original code) */
    /* For production, could implement Fermi distribution with skin */
    double u = random_uniform();
    return ws->R * pow(u, 1.0/3.0);

    /* Alternative: Woods-Saxon rejection sampling
    while (1) {
        double r_max = ws->R + 3.0 * ws->skin_thickness;
        double r = r_max * pow(random_uniform(), 1.0/3.0);
        double rho_r = woods_saxon_density(r, ws);
        double rho_max = ws->rho0;
        if (random_uniform() < pow(rho_r/rho_max, 3.0)) {
            return r;
        }
    }
    */
}

/*
 * Initialize nuclear positions
 */
int initialize_nuclear_positions(
    Particle *particles,
    const NucleusInitConfig *config,
    int start_idx
) {
    /* Calculate minimum separation (from original code) */
    double closest_approach = pow(3.0 / (4.0 * M_PI * config->ws.rho0), 1.0/3.0);
    closest_approach /= sqrt(2.0);  /* Use RMS value */
    double closest_approach_sq = square(closest_approach);

    int particles_created = 0;

    /* Loop over ensembles */
    for (int ens = 0; ens < config->N_ensembles; ens++) {
        /* Loop over nucleons in this nucleus */
        for (int nuc = 0; nuc < config->A; nuc++) {
            int idx = start_idx + particles_created;
            Particle *p = &particles[idx];

            /* Determine if neutron or proton */
            /* First Z particles are protons, rest are neutrons */
            int is_proton = (nuc < config->Z);

            /* Set local radius (same for n and p in simple case) */
            double local_radius = config->ws.R;

            /* Attempt to place particle */
            int too_close = 1;
            int attempts = 0;

            while (too_close && attempts < MAX_INIT_ATTEMPTS) {
                too_close = 0;
                attempts++;

                /* Sample radius and angles */
                double r = sample_woods_saxon_radius(&config->ws);
                double cos_theta = 2.0 * random_uniform() - 1.0;
                double phi = 2.0 * M_PI * random_uniform();
                double sin_theta = sqrt(1.0 - square(cos_theta));

                /* Set position */
                double x_tmp = r * sin_theta * cos(phi);
                double y_tmp = r * sin_theta * sin(phi);
                double z_tmp = r * cos_theta;

                /* Check separation from other nucleons in same ensemble */
                int start_check = start_idx + ens * config->A;
                for (int i = start_check; i < idx; i++) {
                    double dx = x_tmp - particles[i].x;
                    double dy = y_tmp - particles[i].y;
                    double dz = z_tmp - particles[i].z;
                    double r_sq = square(dx) + square(dy) + square(dz);

                    if (r_sq < closest_approach_sq) {
                        too_close = 1;
                        break;
                    }
                }

                if (!too_close) {
                    /* Accept this position */
                    p->x = x_tmp;
                    p->y = y_tmp;
                    p->z = z_tmp;
                }
            }

            if (attempts >= MAX_INIT_ATTEMPTS) {
                fprintf(stderr, "Warning: Could not find valid position after %d attempts\n",
                        MAX_INIT_ATTEMPTS);
                /* Use last attempted position anyway */
            }

            /* Initialize other particle properties */
            p->charge = is_proton ? 1 : 0;
            p->ensemble_id = ens;
            p->nucleon_id = nuc;
            p->is_projectile = config->is_projectile;
            p->num_collisions = 0;

            /* Momenta initialized to zero (set later) */
            p->px = 0.0;
            p->py = 0.0;
            p->pz = 0.0;

            particles_created++;
        }
    }

    return particles_created;
}

/*
 * Calculate local Fermi momentum with LTF corrections
 * Implements Lenk-Pandharipande prescription
 */
double calculate_local_fermi_momentum(
    double rho_local,
    double rho_nabla2,
    double rho_grad2,
    int is_charged,
    double D_pot
) {
    /* Check for very low density */
    if (rho_local < RHO_MIN_THRESHOLD) {
        return 0.0;
    }

    /* Calculate correction terms */
    double term1 = 6.0 * rho_local - rho_nabla2;
    double term2 = rho_grad2;

    /* Handle very small term1 (can be slightly negative) */
    if (term1 < RHO_MIN_THRESHOLD) {
        term1 = 0.0;
    } else {
        term1 *= 3.0 / pow(rho_local, 1.0/3.0);
    }

    /* Scale term2 */
    term2 /= pow(rho_local, 4.0/3.0);

    /* Calculate p_F^2 with LTF corrections */
    double pf_sq;

    if (D_pot != 0.0) {
        /* Isospin-dependent case */
        pf_sq = pow(rho_local * 6.0 * M_PI * M_PI / 2.0, 2.0/3.0);
    } else {
        /* Isospin-independent case */
        pf_sq = pow(rho_local * 6.0 * M_PI * M_PI / 4.0, 2.0/3.0);
    }

    /* Apply Lenk-Pandharipande correction */
    double correction = pow(6.0 * M_PI * M_PI / 2.0, 2.0/3.0) * (term1 - term2) * A_FF;
    pf_sq -= correction;

    /* Check for negative p_F^2 (can happen at very low density) */
    if (pf_sq < 0.0) {
        /* Fall back to uncorrected value */
        pf_sq = pow(rho_local * 6.0 * M_PI * M_PI / 2.0, 2.0/3.0);
        if (pf_sq < 0.0) {
            return 0.0;
        }
    }

    /* Return p_F in MeV (HBARC is in GeV·fm, so multiply by 1000) */
    return sqrt(pf_sq) * HBARC * 1000.0;
}

/*
 * Initialize nuclear momenta using LTF
 */
int initialize_nuclear_momenta(
    Particle *particles,
    int num_particles,
    const Grid *grid,
    double D_pot
) {
    /* Grid parameters */
    double dx_inv = 1.0 / grid->dx;

    for (int i = 0; i < num_particles; i++) {
        Particle *p = &particles[i];

        /* Find grid location */
        int n = where_am_I(p->x, grid->origin_x, dx_inv);
        int m = where_am_I(p->y, grid->origin_y, dx_inv);
        int l = where_am_I(p->z, grid->origin_z, dx_inv);

        /* Check bounds */
        if (n < 2 || n >= grid->nx - 2 ||
            m < 2 || m >= grid->ny - 2 ||
            l < 2 || l >= grid->nz - 2) {
            fprintf(stderr, "Warning: Particle %d outside grid bounds\n", i);
            continue;
        }

        /* Get local density */
        int idx = n + grid->nx * (m + grid->ny * l);
        double rho_local;

        if (D_pot != 0.0) {
            /* Isospin-dependent: separate neutron/proton densities */
            if (p->charge) {
                rho_local = grid->charge_density[idx];
            } else {
                rho_local = grid->baryon_density[idx] - grid->charge_density[idx];
            }
        } else {
            /* Isospin-independent */
            rho_local = grid->baryon_density[idx];
        }

        /* Calculate density derivatives for LTF */
        /* Laplacian (nabla^2 rho) using finite differences */
        double rho_nabla2;
        if (p->charge && D_pot != 0.0) {
            /* Proton */
            int idx_xp = (n+2) + grid->nx * (m + grid->ny * l);
            int idx_xm = (n-2) + grid->nx * (m + grid->ny * l);
            int idx_yp = n + grid->nx * ((m+2) + grid->ny * l);
            int idx_ym = n + grid->nx * ((m-2) + grid->ny * l);
            int idx_zp = n + grid->nx * (m + grid->ny * (l+2));
            int idx_zm = n + grid->nx * (m + grid->ny * (l-2));

            rho_nabla2 = grid->charge_density[idx_xp] + grid->charge_density[idx_xm]
                       + grid->charge_density[idx_yp] + grid->charge_density[idx_ym]
                       + grid->charge_density[idx_zp] + grid->charge_density[idx_zm];
        } else if (!p->charge && D_pot != 0.0) {
            /* Neutron */
            int idx_xp = (n+2) + grid->nx * (m + grid->ny * l);
            int idx_xm = (n-2) + grid->nx * (m + grid->ny * l);
            int idx_yp = n + grid->nx * ((m+2) + grid->ny * l);
            int idx_ym = n + grid->nx * ((m-2) + grid->ny * l);
            int idx_zp = n + grid->nx * (m + grid->ny * (l+2));
            int idx_zm = n + grid->nx * (m + grid->ny * (l-2));

            rho_nabla2 = (grid->baryon_density[idx_xp] - grid->charge_density[idx_xp])
                       + (grid->baryon_density[idx_xm] - grid->charge_density[idx_xm])
                       + (grid->baryon_density[idx_yp] - grid->charge_density[idx_yp])
                       + (grid->baryon_density[idx_ym] - grid->charge_density[idx_ym])
                       + (grid->baryon_density[idx_zp] - grid->charge_density[idx_zp])
                       + (grid->baryon_density[idx_zm] - grid->charge_density[idx_zm]);
        } else {
            /* Isospin-independent */
            int idx_xp = (n+2) + grid->nx * (m + grid->ny * l);
            int idx_xm = (n-2) + grid->nx * (m + grid->ny * l);
            int idx_yp = n + grid->nx * ((m+2) + grid->ny * l);
            int idx_ym = n + grid->nx * ((m-2) + grid->ny * l);
            int idx_zp = n + grid->nx * (m + grid->ny * (l+2));
            int idx_zm = n + grid->nx * (m + grid->ny * (l-2));

            rho_nabla2 = grid->baryon_density[idx_xp] + grid->baryon_density[idx_xm]
                       + grid->baryon_density[idx_yp] + grid->baryon_density[idx_ym]
                       + grid->baryon_density[idx_zp] + grid->baryon_density[idx_zm];
        }

        /* Gradient squared */
        double rho_grad2;
        if (p->charge && D_pot != 0.0) {
            int idx_xp = (n+1) + grid->nx * (m + grid->ny * l);
            int idx_xm = (n-1) + grid->nx * (m + grid->ny * l);
            int idx_yp = n + grid->nx * ((m+1) + grid->ny * l);
            int idx_ym = n + grid->nx * ((m-1) + grid->ny * l);
            int idx_zp = n + grid->nx * (m + grid->ny * (l+1));
            int idx_zm = n + grid->nx * (m + grid->ny * (l-1));

            double dx_rho = grid->charge_density[idx_xp] - grid->charge_density[idx_xm];
            double dy_rho = grid->charge_density[idx_yp] - grid->charge_density[idx_ym];
            double dz_rho = grid->charge_density[idx_zp] - grid->charge_density[idx_zm];
            rho_grad2 = square(dx_rho) + square(dy_rho) + square(dz_rho);
        } else if (!p->charge && D_pot != 0.0) {
            int idx_xp = (n+1) + grid->nx * (m + grid->ny * l);
            int idx_xm = (n-1) + grid->nx * (m + grid->ny * l);
            int idx_yp = n + grid->nx * ((m+1) + grid->ny * l);
            int idx_ym = n + grid->nx * ((m-1) + grid->ny * l);
            int idx_zp = n + grid->nx * (m + grid->ny * (l+1));
            int idx_zm = n + grid->nx * (m + grid->ny * (l-1));

            double dx_rho = (grid->baryon_density[idx_xp] - grid->charge_density[idx_xp])
                          - (grid->baryon_density[idx_xm] - grid->charge_density[idx_xm]);
            double dy_rho = (grid->baryon_density[idx_yp] - grid->charge_density[idx_yp])
                          - (grid->baryon_density[idx_ym] - grid->charge_density[idx_ym]);
            double dz_rho = (grid->baryon_density[idx_zp] - grid->charge_density[idx_zp])
                          - (grid->baryon_density[idx_zm] - grid->charge_density[idx_zm]);
            rho_grad2 = square(dx_rho) + square(dy_rho) + square(dz_rho);
        } else {
            int idx_xp = (n+1) + grid->nx * (m + grid->ny * l);
            int idx_xm = (n-1) + grid->nx * (m + grid->ny * l);
            int idx_yp = n + grid->nx * ((m+1) + grid->ny * l);
            int idx_ym = n + grid->nx * ((m-1) + grid->ny * l);
            int idx_zp = n + grid->nx * (m + grid->ny * (l+1));
            int idx_zm = n + grid->nx * (m + grid->ny * (l-1));

            double dx_rho = grid->baryon_density[idx_xp] - grid->baryon_density[idx_xm];
            double dy_rho = grid->baryon_density[idx_yp] - grid->baryon_density[idx_ym];
            double dz_rho = grid->baryon_density[idx_zp] - grid->baryon_density[idx_zm];
            rho_grad2 = square(dx_rho) + square(dy_rho) + square(dz_rho);
        }

        /* Calculate local Fermi momentum */
        double p_fermi = calculate_local_fermi_momentum(
            rho_local, rho_nabla2, rho_grad2, p->charge, D_pot
        );

        /* Sample momentum uniformly in Fermi sphere */
        double u = pow(random_uniform(), 1.0/3.0);
        double p_mag = p_fermi * u;

        /* Random angles */
        double cos_theta = 2.0 * random_uniform() - 1.0;
        double phi = 2.0 * M_PI * random_uniform();
        double sin_theta = sqrt(1.0 - square(cos_theta));

        /* Set momentum (in MeV) */
        p->px = p_mag * sin_theta * cos(phi);
        p->py = p_mag * sin_theta * sin(phi);
        p->pz = p_mag * cos_theta;
    }

    return 0;
}

/*
 * Zero center-of-mass momentum
 */
void zero_cm_momentum(
    Particle *particles,
    int start_idx,
    int count,
    int N_ensembles
) {
    /* Calculate average momentum */
    double px_avg = 0.0, py_avg = 0.0, pz_avg = 0.0;

    for (int i = start_idx; i < start_idx + count; i++) {
        px_avg += particles[i].px;
        py_avg += particles[i].py;
        pz_avg += particles[i].pz;
    }

    double norm = 1.0 / count;
    px_avg *= norm;
    py_avg *= norm;
    pz_avg *= norm;

    /* Shift all momenta */
    for (int i = start_idx; i < start_idx + count; i++) {
        particles[i].px -= px_avg;
        particles[i].py -= py_avg;
        particles[i].pz -= pz_avg;
    }
}

/*
 * Generate full nucleus
 */
int generate_nucleus(
    Particle *particles,
    const NucleusInitConfig *config,
    const Grid *grid,
    int start_idx,
    double D_pot
) {
    /* Initialize positions */
    int n_created = initialize_nuclear_positions(particles, config, start_idx);

    if (n_created != config->A * config->N_ensembles) {
        fprintf(stderr, "Error: Created %d particles, expected %d\n",
                n_created, config->A * config->N_ensembles);
        return -1;
    }

    /* Note: Momenta should be initialized after density grids are computed */
    /* This is typically done in a separate step in the full initialization */

    return 0;
}
