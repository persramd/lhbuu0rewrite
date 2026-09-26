#include "performance.h"
#include "utils.h"
#include <stdio.h>

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#ifdef __linux__
#include <unistd.h>
#endif

size_t perf_get_l2_cache_size(void) {
#ifdef __APPLE__
    size_t cache_size = 0;
    size_t len = sizeof(cache_size);
    if (sysctlbyname("hw.l2cachesize", &cache_size, &len, NULL, 0) == 0) {
        return cache_size;
    }
#endif
#ifdef __linux__
    FILE *fp = fopen("/sys/devices/system/cpu/cpu0/cache/index2/size", "r");
    if (fp) {
        int size_kb;
        if (fscanf(fp, "%dK", &size_kb) == 1) {
            fclose(fp);
            return size_kb * 1024;
        }
        fclose(fp);
    }
#endif
    /* Default estimate: 512KB per core, 4 cores sharing */
    return 2 * 1024 * 1024;  /* 2 MB */
}

size_t perf_get_l3_cache_size(void) {
#ifdef __APPLE__
    size_t cache_size = 0;
    size_t len = sizeof(cache_size);
    if (sysctlbyname("hw.l3cachesize", &cache_size, &len, NULL, 0) == 0) {
        return cache_size;
    }
    return 0;  /* M-series doesn't have L3, just large L2 */
#endif
#ifdef __linux__
    FILE *fp = fopen("/sys/devices/system/cpu/cpu0/cache/index3/size", "r");
    if (fp) {
        int size_kb;
        if (fscanf(fp, "%dK", &size_kb) == 1) {
            fclose(fp);
            return size_kb * 1024;
        }
        fclose(fp);
    }
#endif
    return 8 * 1024 * 1024;  /* 8 MB default */
}

int perf_get_num_cores(void) {
#ifdef __APPLE__
    int num_cores = 0;
    size_t len = sizeof(num_cores);
    if (sysctlbyname("hw.physicalcpu", &num_cores, &len, NULL, 0) == 0) {
        return num_cores;
    }
#endif
#ifdef __linux__
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
    return 4;  /* Default */
}

PerformanceMode perf_analyze_memory(const SimConfig *config, MemoryProfile *profile) {
    /* Calculate memory requirements */
    size_t n_particles = config->num_ensembles *
                        (config->projectile.A + config->target.A);

    size_t grid_cells = config->grid.nx * config->grid.ny * config->grid.nz;

    /* Particle data */
    profile->particle_memory = n_particles * sizeof(Particle);

    /* Grid data (4 float grids + cell tracking) */
    profile->grid_memory = grid_cells * 4 * sizeof(float) +  /* Density/potential grids */
                          grid_cells * sizeof(int);          /* Cell counts */

    /* Neighbor lists (estimate 100 neighbors per particle) */
    profile->neighbor_memory = n_particles * 100 * (sizeof(int) + sizeof(float));

    /* Total working set */
    profile->total_memory = profile->particle_memory +
                           profile->grid_memory +
                           profile->neighbor_memory;

    /* Get cache sizes */
    profile->l2_cache_size = perf_get_l2_cache_size();
    profile->l3_cache_size = perf_get_l3_cache_size();

    /* Determine mode */
    if (profile->total_memory < profile->l2_cache_size * 0.8) {
        profile->mode = PERF_CACHE_FIT;
    } else if (profile->l3_cache_size > 0 &&
               profile->total_memory < profile->l3_cache_size * 0.8) {
        profile->mode = PERF_CACHE_AWARE;
    } else {
        profile->mode = PERF_MEMORY_BOUND;
    }

    return profile->mode;
}

void perf_print_profile(const MemoryProfile *profile) {
    printf("\n=== Memory and Performance Profile ===\n");
    printf("Working set size:\n");
    printf("  Particles:  %6.2f MB\n", profile->particle_memory / (1024.0 * 1024.0));
    printf("  Grid:       %6.2f MB\n", profile->grid_memory / (1024.0 * 1024.0));
    printf("  Neighbors:  %6.2f MB\n", profile->neighbor_memory / (1024.0 * 1024.0));
    printf("  Total:      %6.2f MB\n", profile->total_memory / (1024.0 * 1024.0));
    printf("\nCache hierarchy:\n");
    printf("  L2 cache:   %6.2f MB\n", profile->l2_cache_size / (1024.0 * 1024.0));
    if (profile->l3_cache_size > 0) {
        printf("  L3 cache:   %6.2f MB\n", profile->l3_cache_size / (1024.0 * 1024.0));
    }
    printf("  CPU cores:  %d\n", perf_get_num_cores());

    printf("\nPerformance mode: ");
    switch (profile->mode) {
        case PERF_CACHE_FIT:
            printf("CACHE_FIT (everything in L2)\n");
            printf("  Strategy: Simple algorithms, minimal overhead\n");
            break;
        case PERF_CACHE_AWARE:
            printf("CACHE_AWARE (fits in L3)\n");
            printf("  Strategy: Periodic sorting, neighbor lists\n");
            break;
        case PERF_MEMORY_BOUND:
            printf("MEMORY_BOUND (large system)\n");
            printf("  Strategy: Aggressive optimization, frequent sorting\n");
            break;
    }
    printf("======================================\n\n");
}
