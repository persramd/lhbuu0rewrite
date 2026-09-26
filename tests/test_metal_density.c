#include "grid.h"
#include "gpu_compute_metal.h"
#include "types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* Compare two density grids */
static int compare_grids(const float *grid1, const float *grid2, size_t n_cells,
                         float tolerance, const char *name) {
    float max_diff = 0.0f;
    float avg_diff = 0.0f;
    int mismatch_count = 0;

    for (size_t i = 0; i < n_cells; i++) {
        float diff = fabsf(grid1[i] - grid2[i]);
        if (diff > tolerance) {
            mismatch_count++;
        }
        if (diff > max_diff) max_diff = diff;
        avg_diff += diff;
    }

    avg_diff /= n_cells;

    printf("\n%s comparison:\n", name);
    printf("  Max difference: %e\n", max_diff);
    printf("  Avg difference: %e\n", avg_diff);
    printf("  Mismatches (>%.0e): %d / %zu\n", tolerance, mismatch_count, n_cells);

    return mismatch_count == 0;
}

int main(void) {
    printf("=== Metal Grid Density Benchmark ===\n\n");

    /* Grid parameters (64^3 cells) */
    int nx = 64, ny = 64, nz = 64;
    float dx = 0.5f;  /* fm */
    int ff_range = 2;

    /* Number of particles (typical for Au+Au collision) */
    size_t n_particles = 400;

    printf("Configuration:\n");
    printf("  Grid: %dx%dx%d = %d cells\n", nx, ny, nz, nx*ny*nz);
    printf("  Cell size: %.2f fm\n", dx);
    printf("  Particles: %zu\n", n_particles);
    printf("  Form factor range: ±%d cells\n\n", ff_range);

    /* Initialize grids */
    Grid grid_cpu, grid_gpu;
    if (!grid_init(&grid_cpu, nx, ny, nz, dx)) {
        fprintf(stderr, "Failed to initialize CPU grid\n");
        return 1;
    }
    if (!grid_init(&grid_gpu, nx, ny, nz, dx)) {
        fprintf(stderr, "Failed to initialize GPU grid\n");
        return 1;
    }

    /* Generate random particles */
    Particle *particles = malloc(n_particles * sizeof(Particle));
    if (!particles) {
        fprintf(stderr, "Failed to allocate particles\n");
        return 1;
    }

    srand(42);  /* Reproducible random seed */
    float box_size = nx * dx;
    for (size_t i = 0; i < n_particles; i++) {
        particles[i].x = (rand() / (float)RAND_MAX - 0.5f) * box_size * 0.8f;
        particles[i].y = (rand() / (float)RAND_MAX - 0.5f) * box_size * 0.8f;
        particles[i].z = (rand() / (float)RAND_MAX - 0.5f) * box_size * 0.8f;
        particles[i].px = 0.0f;
        particles[i].py = 0.0f;
        particles[i].pz = 0.0f;
        particles[i].charge = rand() % 2;  /* 50% protons */
        particles[i].ensemble_id = 0;
        particles[i].nucleon_id = i;
        particles[i].is_projectile = 0;
        particles[i].num_collisions = 0;
        particles[i].grid_cell = -1;
    }

    /* Initialize Metal */
    printf("Initializing Metal...\n");
    if (!gpu_init_metal()) {
        fprintf(stderr, "Failed to initialize Metal\n");
        return 1;
    }

    /* Warm-up run (compile kernels, etc.) */
    printf("\nWarm-up run...\n");
    grid_calculate_density_metal(&grid_gpu, particles, n_particles, ff_range);

    /* Benchmark CPU version */
    printf("\nBenchmarking CPU version...\n");
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    const int n_iterations = 10;
    for (int i = 0; i < n_iterations; i++) {
        grid_calculate_density(&grid_cpu, particles, n_particles, ff_range);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double cpu_time = (end.tv_sec - start.tv_sec) +
                      (end.tv_nsec - start.tv_nsec) / 1e9;
    cpu_time /= n_iterations;

    printf("CPU time: %.3f ms\n", cpu_time * 1000.0);

    /* Benchmark GPU version */
    printf("\nBenchmarking GPU version...\n");
    clock_gettime(CLOCK_MONOTONIC, &start);

    for (int i = 0; i < n_iterations; i++) {
        grid_calculate_density_metal(&grid_gpu, particles, n_particles, ff_range);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double gpu_time = (end.tv_sec - start.tv_sec) +
                      (end.tv_nsec - start.tv_nsec) / 1e9;
    gpu_time /= n_iterations;

    printf("GPU time: %.3f ms\n", gpu_time * 1000.0);

    /* Performance summary */
    printf("\n=== Performance Summary ===\n");
    printf("CPU time:    %.3f ms\n", cpu_time * 1000.0);
    printf("GPU time:    %.3f ms\n", gpu_time * 1000.0);
    printf("Speedup:     %.2fx\n", cpu_time / gpu_time);
    printf("Throughput:  %.1f particles/ms (GPU)\n", n_particles / (gpu_time * 1000.0));

    /* Verify correctness */
    size_t n_cells = nx * ny * nz;
    int baryon_match = compare_grids(grid_cpu.baryon_density,
                                     grid_gpu.baryon_density,
                                     n_cells, 1e-5f, "Baryon density");

    int charge_match = compare_grids(grid_cpu.charge_density,
                                     grid_gpu.charge_density,
                                     n_cells, 1e-5f, "Charge density");

    /* Sample output */
    printf("\nSample density values (cell [32,32,32]):\n");
    int idx = grid_index(&grid_cpu, 32, 32, 32);
    printf("  CPU baryon: %e\n", grid_cpu.baryon_density[idx]);
    printf("  GPU baryon: %e\n", grid_gpu.baryon_density[idx]);
    printf("  CPU charge: %e\n", grid_cpu.charge_density[idx]);
    printf("  GPU charge: %e\n", grid_gpu.charge_density[idx]);

    /* Cleanup */
    gpu_cleanup_metal();
    grid_free(&grid_cpu);
    grid_free(&grid_gpu);
    free(particles);

    if (baryon_match && charge_match) {
        printf("\n*** PASS: CPU and GPU results match ***\n");
        return 0;
    } else {
        printf("\n*** FAIL: CPU and GPU results differ ***\n");
        return 1;
    }
}
