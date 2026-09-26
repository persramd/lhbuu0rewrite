#include "config.h"
#include "utils.h"
#include <string.h>
#include <stdlib.h>

void config_set_defaults(SimConfig *config) {
    memset(config, 0, sizeof(SimConfig));

    /* Projectile: proton */
    config->projectile.Z = 1;
    config->projectile.A = 1;

    /* Target: Ca-40 */
    config->target.Z = 20;
    config->target.A = 40;

    /* Kinematics */
    config->beam_energy = 200.0f;  /* MeV/nucleon */
    config->impact_parameter = 0.5f;  /* Normalized */

    /* Time evolution */
    config->dt = 0.2f;  /* fm/c */
    config->t_max = 100.0f;  /* fm/c */
    config->diagnostic_interval = 50;  /* steps */

    /* Ensembles */
    config->num_ensembles = 100;

    /* Grid */
    config->grid.nx = 60;
    config->grid.ny = 60;
    config->grid.nz = 60;
    config->grid.dx = 1.0f;  /* fm */
    config->form_factor_range = 2;
    config->sort_interval = 10;

    /* Potential: S_K200 (soft, compressibility=200 MeV) */
    config->potential.A = -0.356f;  /* GeV */
    config->potential.B = 0.303f;   /* GeV */
    config->potential.sigma = 1.167f;
    config->potential.C = 0.0f;  /* No momentum-dependence by default */
    config->potential.lambda = 0.0f;
    config->potential.D = 0.032f;  /* Isospin */
    config->potential.A_surf = 0.0f;
    config->potential.lambda_surf = 0.0f;
    config->potential.rho0 = 0.16f;  /* fm^-3 */
    config->potential.enable_coulomb = 1;
    config->potential.enable_momentum_dep = 0;

    /* Collisions */
    config->collision.sigma_nn_max = 40.0f;  /* mb */
    config->collision.sigma_rho_reduction = 0.2f;
    config->collision.enable_pauli_blocking = 1;
    config->collision.pauli_ensemble_count = 100;

    /* Output */
    strcpy(config->output_prefix, "lhbuu");
    strcpy(config->output_path, "./");
    config->full_output = 1;
    config->animation_interval = 0;  /* Disabled */
}

int config_load(const char *filename, SimConfig *config) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        warning("Cannot open config file '%s', using defaults", filename);
        config_set_defaults(config);
        return 0;
    }

    /* Start with defaults */
    config_set_defaults(config);

    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        /* Skip comments and empty lines */
        if (line[0] == '#' || line[0] == '\n') continue;

        char key[64], value[192];
        if (sscanf(line, "%63s = %191[^\n]", key, value) == 2) {
            /* Parse configuration */
            if (strcmp(key, "projectile_Z") == 0) config->projectile.Z = atoi(value);
            else if (strcmp(key, "projectile_A") == 0) config->projectile.A = atoi(value);
            else if (strcmp(key, "target_Z") == 0) config->target.Z = atoi(value);
            else if (strcmp(key, "target_A") == 0) config->target.A = atoi(value);
            else if (strcmp(key, "beam_energy") == 0) config->beam_energy = atof(value);
            else if (strcmp(key, "impact_parameter") == 0) config->impact_parameter = atof(value);
            else if (strcmp(key, "dt") == 0) config->dt = atof(value);
            else if (strcmp(key, "t_max") == 0) config->t_max = atof(value);
            else if (strcmp(key, "num_ensembles") == 0) config->num_ensembles = atoi(value);
            else if (strcmp(key, "output_prefix") == 0) strncpy(config->output_prefix, value, 255);
            else if (strcmp(key, "potential_A") == 0) config->potential.A = atof(value);
            else if (strcmp(key, "potential_B") == 0) config->potential.B = atof(value);
            else if (strcmp(key, "potential_sigma") == 0) config->potential.sigma = atof(value);
        }
    }

    fclose(fp);
    return 1;
}

int config_validate(const SimConfig *config) {
    if (config->projectile.A <= 0 || config->target.A <= 0) {
        error_exit("Invalid nuclear mass numbers");
    }
    if (config->beam_energy <= 0) {
        error_exit("Beam energy must be positive");
    }
    if (config->dt <= 0 || config->t_max <= 0) {
        error_exit("Time parameters must be positive");
    }
    if (config->grid.nx <= 0 || config->grid.ny <= 0 || config->grid.nz <= 0) {
        error_exit("Grid dimensions must be positive");
    }
    return 1;
}

void config_print(FILE *fp, const SimConfig *config) {
    fprintf(fp, "# LHBUU Configuration\n");
    fprintf(fp, "projectile_Z = %d\n", config->projectile.Z);
    fprintf(fp, "projectile_A = %d\n", config->projectile.A);
    fprintf(fp, "target_Z = %d\n", config->target.Z);
    fprintf(fp, "target_A = %d\n", config->target.A);
    fprintf(fp, "beam_energy = %.2f  # MeV/nucleon\n", config->beam_energy);
    fprintf(fp, "impact_parameter = %.3f\n", config->impact_parameter);
    fprintf(fp, "dt = %.2f  # fm/c\n", config->dt);
    fprintf(fp, "t_max = %.2f  # fm/c\n", config->t_max);
    fprintf(fp, "num_ensembles = %d\n", config->num_ensembles);
    fprintf(fp, "potential_A = %.4f  # GeV\n", config->potential.A);
    fprintf(fp, "potential_B = %.4f  # GeV\n", config->potential.B);
    fprintf(fp, "potential_sigma = %.4f\n", config->potential.sigma);
}
