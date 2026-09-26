#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "types.h"
#include "config.h"
#include "initialize.h"
#include "integrator.h"
#include "collision.h"
#include "diagnostics.h"
#include "performance.h"
#include "utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void print_usage(const char *progname) {
    printf("Usage: %s [options]\n", progname);
    printf("Options:\n");
    printf("  -c <file>     Configuration file\n");
    printf("  -e <energy>   Beam energy (MeV/nucleon)\n");
    printf("  -b <impact>   Impact parameter (0-1)\n");
    printf("  -t <time>     Max time (fm/c)\n");
    printf("  -o <prefix>   Output prefix\n");
    printf("  -h            Show this help\n");
}

int main(int argc, char *argv[]) {
    printf("LHBUU - Lattice Hamiltonian BUU Transport Code (Rewrite)\n\n");

    SimConfig config;
    SimState state;
    char config_file[256] = "";
    int use_config_file = 0;

    /* Parse command line */
    config_set_defaults(&config);

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-c") == 0 && i+1 < argc) {
            strncpy(config_file, argv[++i], 255);
            use_config_file = 1;
        } else if (strcmp(argv[i], "-e") == 0 && i+1 < argc) {
            config.beam_energy = atof(argv[++i]);
        } else if (strcmp(argv[i], "-b") == 0 && i+1 < argc) {
            config.impact_parameter = atof(argv[++i]);
        } else if (strcmp(argv[i], "-t") == 0 && i+1 < argc) {
            config.t_max = atof(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0 && i+1 < argc) {
            strncpy(config.output_prefix, argv[++i], 255);
        }
    }

    /* Load config file if specified */
    if (use_config_file) {
        config_load(config_file, &config);
    }

    /* Validate configuration */
    if (!config_validate(&config)) {
        return 1;
    }

    /* Print configuration */
    printf("Configuration:\n");
    printf("  Projectile: %d nucleons (Z=%d)\n", config.projectile.A, config.projectile.Z);
    printf("  Target: %d nucleons (Z=%d)\n", config.target.A, config.target.Z);
    printf("  Beam energy: %.2f MeV/nucleon\n", config.beam_energy);
    printf("  Impact parameter: %.3f\n", config.impact_parameter);
    printf("  Time: dt=%.2f fm/c, t_max=%.2f fm/c\n", config.dt, config.t_max);
    printf("  Grid: %dx%dx%d, dx=%.2f fm\n", config.grid.nx, config.grid.ny, config.grid.nz, config.grid.dx);
    printf("  Ensembles: %d\n", config.num_ensembles);
    printf("\n");

    /* Analyze memory and performance */
    MemoryProfile mem_profile;
    perf_analyze_memory(&config, &mem_profile);
    perf_print_profile(&mem_profile);

    /* Set OpenMP threads if available */
    #ifdef _OPENMP
    int num_threads = perf_get_num_cores();
    omp_set_num_threads(num_threads);
    printf("OpenMP enabled: %d threads\n\n", num_threads);
    #endif

    /* Initialize random number generator */
    random_init((unsigned int)time(NULL));

    /* Initialize simulation */
    printf("Initializing simulation...\n");
    if (!sim_init(&state, &config)) {
        error_exit("Simulation initialization failed");
    }
    printf("  Total particles: %zu\n", state.num_particles);
    printf("\n");

    /* Open diagnostics file */
    char diag_filename[512];
    snprintf(diag_filename, sizeof(diag_filename), "%s%s_diag.txt",
             config.output_path, config.output_prefix);
    FILE *diag_fp = file_open(diag_filename, "w");
    fprintf(diag_fp, "# timestep time KE PE E px py pz collisions\n");

    /* Initial diagnostics */
    diagnostics_write_timestep(diag_fp, &state);

    /* Main simulation loop */
    printf("Running simulation...\n");
    int max_steps = (int)(config.t_max / config.dt);

    for (int step = 0; step < max_steps; step++) {
        /* Time integration */
        integrator_step(&state);

        /* Collisions */
        if (config.collision.sigma_nn_max > 0) {
            collision_process_all(&state);
        }

        /* Diagnostics */
        if (step % config.diagnostic_interval == 0) {
            diagnostics_write_timestep(diag_fp, &state);
            printf("  Step %d/%d (t=%.2f fm/c), E=%.4f GeV, collisions=%d\n",
                   step, max_steps, state.current_time,
                   diagnostics_total_energy(&state),
                   state.total_collisions);
        }
    }

    printf("Simulation complete.\n\n");

    /* Final diagnostics */
    diagnostics_print_summary(&state);
    fclose(diag_fp);

    /* Write final phase space */
    if (config.full_output) {
        char ps_filename[512];
        snprintf(ps_filename, sizeof(ps_filename), "%s%s_final.dat",
                 config.output_path, config.output_prefix);
        diagnostics_write_phase_space(ps_filename, &state);
        printf("Phase space written to: %s\n", ps_filename);
    }

    printf("Diagnostics written to: %s\n", diag_filename);

    /* Cleanup */
    sim_free(&state);

    return 0;
}
