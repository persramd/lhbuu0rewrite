#include "gpu_compute.h"
#include "grid.h"
#include "types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* Compare two density arrays */
static void compare_results(const float *cpu_data, const float *gpu_data,
                            size_t n, const char *name) {
    float max_diff = 0.0f;
    float avg_diff = 0.0f;
    int mismatch_count = 0;

    for (size_t i = 0; i < n; i++) {
        float diff = fabsf(cpu_data[i] - gpu_data[i]);
        if (diff > 1e-5f) {
            mismatch_count++;
        }
        if (diff > max_diff) max_diff = diff;
        avg_diff += diff;
    }

    avg_diff /= n;

    printf("\n%s comparison:\n", name);
    printf("  Max difference: %e\n", max_diff);
    printf("  Avg difference: %e\n", avg_diff);
    printf("  Mismatches (>1e-5): %d / %zu (%.2f%%)\n",
           mismatch_count, n, 100.0f * mismatch_count / n);
}

int main(void) {
    printf("=== LHBUU Metal GPU Density Calculation Test ===\n\n");

    /* Grid configuration (64³ cells as specified) */
    int nx = 64, ny = 64, nz = 64;
    float dx = 0.5f;  /* fm */
    int ff_range = 2;  /* ±2 cells form factor support */

    /* Particle count (typical for Au+Au) */
    size_t n_particles = 400;

    printf("Configuration:\n");
    printf("  Grid: %dx%dx%d = %d cells\n", nx, ny, nz, nx*ny*nz);
    printf("  Cell size: %.2f fm\n", dx);
    printf("  Particles: %zu\n", n_particles);
    printf("  Form factor range: ±%d cells\n\n", ff_range);

    /* Initialize GPU */
    printf("Initializing Metal GPU...\n");
    GPUContext *gpu_ctx = gpu_init();

    if (!gpu_ctx || !gpu_is_available(gpu_ctx)) {
        fprintf(stderr, "Failed to initialize GPU: %s\n",
                gpu_ctx ? gpu_get_error(gpu_ctx) : "NULL context");
        return 1;
    }

    GPUDeviceInfo info = gpu_get_device_info(gpu_ctx);
    printf("  Device: %s\n", info.device_name);
    printf("  Backend: %s\n", info.backend == GPU_BACKEND_METAL ? "Metal" : "Unknown");
    printf("  Memory: %.1f GB %s\n", info.total_memory / 1e9,
           info.unified_memory ? "(unified)" : "");
    printf("  Compute units: %d cores\n", info.compute_units);
    printf("  Max threads/group: %zu\n\n", info.max_threads_per_threadgroup);

    /* Generate random particles */
    Particle *particles = malloc(n_particles * sizeof(Particle));
    if (!particles) {
        fprintf(stderr, "Failed to allocate particles\n");
        return 1;
    }

    srand(42);  /* Reproducible */
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

    /* Allocate density arrays */
    size_t n_cells = nx * ny * nz;
    float *cpu_baryon = calloc(n_cells, sizeof(float));
    float *cpu_charge = calloc(n_cells, sizeof(float));
    float *gpu_baryon = calloc(n_cells, sizeof(float));
    float *gpu_charge = calloc(n_cells, sizeof(float));
    float *density_combined = calloc(n_cells * 2, sizeof(float));  /* baryon + charge */

    if (!cpu_baryon || !cpu_charge || !gpu_baryon || !gpu_charge || !density_combined) {
        fprintf(stderr, "Failed to allocate density arrays\n");
        return 1;
    }

    float origin_x = -0.5f * nx * dx;
    float origin_y = -0.5f * ny * dx;
    float origin_z = -0.5f * nz * dx;

    /* ===== CPU Reference Calculation ===== */
    printf("Running CPU reference calculation...\n");
    struct timespec cpu_start, cpu_end;
    clock_gettime(CLOCK_MONOTONIC, &cpu_start);

    /* Combine baryon and charge for CPU fallback */
    cpu_calculate_density(cpu_baryon, particles, n_particles,
                         nx, ny, nz, dx,
                         origin_x, origin_y, origin_z, ff_range);

    /* Charge density (separate pass for simplicity in test) */
    memset(cpu_charge, 0, n_cells * sizeof(float));
    for (size_t i = 0; i < n_particles; i++) {
        if (!particles[i].charge) continue;

        int ix = (int)((particles[i].x - origin_x) / dx + 0.5f);
        int iy = (int)((particles[i].y - origin_y) / dx + 0.5f);
        int iz = (int)((particles[i].z - origin_z) / dx + 0.5f);

        for (int dx_off = -ff_range; dx_off <= ff_range; dx_off++) {
            int cx = ix + dx_off;
            if (cx < 0 || cx >= nx) continue;

            for (int dy_off = -ff_range; dy_off <= ff_range; dy_off++) {
                int cy = iy + dy_off;
                if (cy < 0 || cy >= ny) continue;

                for (int dz_off = -ff_range; dz_off <= ff_range; dz_off++) {
                    int cz = iz + dz_off;
                    if (cz < 0 || cz >= nz) continue;

                    float dist2 = (dx_off*dx) * (dx_off*dx) +
                                  (dy_off*dx) * (dy_off*dx) +
                                  (dz_off*dx) * (dz_off*dx);
                    float sigma_ff = ff_range * dx / 2.0f;
                    float norm = 1.0f / (sigma_ff * sqrtf(2.0f * M_PI));
                    float weight = norm * expf(-0.5f * dist2 / (sigma_ff * sigma_ff));

                    int idx = cx * ny * nz + cy * nz + cz;
                    cpu_charge[idx] += weight / (dx * dx * dx);
                }
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &cpu_end);
    double cpu_time = (cpu_end.tv_sec - cpu_start.tv_sec) +
                      (cpu_end.tv_nsec - cpu_start.tv_nsec) / 1e9;

    printf("  CPU time: %.3f ms\n", cpu_time * 1000.0);

    /* ===== GPU Calculation ===== */
    printf("\nRunning GPU calculation...\n");

    /* Create GPU buffers */
    GPUBuffer *particle_buf = gpu_buffer_alloc_with_data(gpu_ctx, particles,
                                                         n_particles * sizeof(Particle));
    GPUBuffer *density_buf = gpu_buffer_alloc(gpu_ctx, n_cells * 2 * sizeof(float));

    if (!particle_buf || !density_buf) {
        fprintf(stderr, "Failed to allocate GPU buffers: %s\n",
                gpu_get_error(gpu_ctx));
        return 1;
    }

    /* Zero density buffer */
    memset(density_buf->host_ptr, 0, n_cells * 2 * sizeof(float));

    struct timespec gpu_start, gpu_end;
    clock_gettime(CLOCK_MONOTONIC, &gpu_start);

    /* Execute density kernel */
    bool success = gpu_calculate_density(
        gpu_ctx,
        density_buf,
        particle_buf,
        n_particles,
        nx, ny, nz, dx,
        origin_x, origin_y, origin_z,
        ff_range
    );

    if (!success) {
        fprintf(stderr, "GPU density calculation failed: %s\n",
                gpu_get_error(gpu_ctx));
        return 1;
    }

    gpu_synchronize(gpu_ctx);

    clock_gettime(CLOCK_MONOTONIC, &gpu_end);
    double gpu_time = (gpu_end.tv_sec - gpu_start.tv_sec) +
                      (gpu_end.tv_nsec - gpu_start.tv_nsec) / 1e9;

    printf("  GPU time: %.3f ms\n", gpu_time * 1000.0);

    /* Copy results back */
    float *density_ptr = (float*)density_buf->host_ptr;
    memcpy(gpu_baryon, density_ptr, n_cells * sizeof(float));
    memcpy(gpu_charge, density_ptr + n_cells, n_cells * sizeof(float));

    /* ===== Performance Summary ===== */
    printf("\n=== Performance Summary ===\n");
    printf("CPU time:    %.3f ms\n", cpu_time * 1000.0);
    printf("GPU time:    %.3f ms\n", gpu_time * 1000.0);
    printf("Speedup:     %.2fx\n", cpu_time / gpu_time);
    printf("Throughput:  %.1f Mparticles/s (GPU)\n",
           n_particles / (gpu_time * 1e6));

    /* ===== Correctness Verification ===== */
    compare_results(cpu_baryon, gpu_baryon, n_cells, "Baryon density");
    compare_results(cpu_charge, gpu_charge, n_cells, "Charge density");

    /* Sample output */
    printf("\nSample values (center cell [32,32,32]):\n");
    int center_idx = 32 * ny * nz + 32 * nz + 32;
    printf("  CPU baryon: %e\n", cpu_baryon[center_idx]);
    printf("  GPU baryon: %e\n", gpu_baryon[center_idx]);
    printf("  CPU charge: %e\n", cpu_charge[center_idx]);
    printf("  GPU charge: %e\n", gpu_charge[center_idx]);

    /* Cleanup */
    gpu_buffer_free(particle_buf);
    gpu_buffer_free(density_buf);
    gpu_cleanup(gpu_ctx);

    free(particles);
    free(cpu_baryon);
    free(cpu_charge);
    free(gpu_baryon);
    free(gpu_charge);
    free(density_combined);

    printf("\n=== Test Complete ===\n");
    return 0;
}
