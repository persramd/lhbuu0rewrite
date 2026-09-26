#include "gpu_compute.h"
#include "types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Simple test for GPU compute infrastructure */

static void print_device_info(const GPUDeviceInfo *info) {
    printf("\nGPU Device Information:\n");
    printf("  Backend: ");
    switch (info->backend) {
        case GPU_BACKEND_NONE:  printf("None (CPU fallback)\n"); break;
        case GPU_BACKEND_METAL: printf("Metal\n"); break;
        case GPU_BACKEND_HIP:   printf("HIP/ROCm\n"); break;
    }
    printf("  Device: %s\n", info->device_name);
    printf("  Compute units: %d\n", info->compute_units);
    printf("  Unified memory: %s\n", info->unified_memory ? "Yes" : "No");
    printf("  Max buffer size: %zu MB\n", info->max_buffer_size / (1024*1024));
    printf("  Total memory: %zu MB\n", info->total_memory / (1024*1024));
    printf("  Max threads/threadgroup: %zu\n", info->max_threads_per_threadgroup);
}

static double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

int main(void) {
    printf("LHBUU GPU Compute Infrastructure Test\n");
    printf("======================================\n");

    /* Initialize GPU context */
    GPUContext *ctx = gpu_init();

    if (!gpu_is_available(ctx)) {
        printf("\nGPU not available - testing CPU fallback\n");
        GPUDeviceInfo info = gpu_get_device_info(ctx);
        print_device_info(&info);
    } else {
        printf("\nGPU initialized successfully\n");
        GPUDeviceInfo info = gpu_get_device_info(ctx);
        print_device_info(&info);
    }

    /* Test parameters */
    const int nx = 64, ny = 64, nz = 64;
    const float dx = 0.5f;
    const float origin = -16.0f;
    const int form_factor_range = 2;
    const size_t num_particles = 1000;

    printf("\nTest configuration:\n");
    printf("  Grid: %d x %d x %d cells\n", nx, ny, nz);
    printf("  Cell spacing: %.2f fm\n", dx);
    printf("  Number of particles: %zu\n", num_particles);
    printf("  Form factor range: %d cells\n", form_factor_range);

    /* Allocate and initialize particle data */
    Particle *particles = malloc(num_particles * sizeof(Particle));
    for (size_t i = 0; i < num_particles; i++) {
        /* Random positions within grid */
        particles[i].x = origin + ((float)rand() / RAND_MAX) * (nx * dx);
        particles[i].y = origin + ((float)rand() / RAND_MAX) * (ny * dx);
        particles[i].z = origin + ((float)rand() / RAND_MAX) * (nz * dx);
        particles[i].charge = rand() % 2;  /* 0=neutron, 1=proton */
    }

    /* Allocate density grid */
    size_t grid_size = nx * ny * nz;
    float *density_cpu = calloc(grid_size, sizeof(float));
    float *density_gpu = calloc(grid_size, sizeof(float));

    /* Test CPU fallback */
    printf("\nTesting CPU implementation...\n");
    double t0 = get_time();
    cpu_calculate_density(density_cpu, particles, num_particles,
                         nx, ny, nz, dx, origin, origin, origin,
                         form_factor_range);
    double cpu_time = get_time() - t0;
    printf("  CPU time: %.3f ms\n", cpu_time * 1000.0);

    /* Calculate some statistics */
    float max_density = 0.0f;
    float total_density = 0.0f;
    for (size_t i = 0; i < grid_size; i++) {
        if (density_cpu[i] > max_density) {
            max_density = density_cpu[i];
        }
        total_density += density_cpu[i];
    }
    printf("  Max density: %.6f fm^-3\n", max_density);
    printf("  Total density: %.6f\n", total_density);

    /* Test GPU implementation if available */
    if (gpu_is_available(ctx)) {
        printf("\nTesting GPU implementation...\n");

        /* Allocate GPU buffers */
        GPUBuffer *particle_buffer = gpu_buffer_alloc_with_data(ctx, particles,
                                                                num_particles * sizeof(Particle));
        GPUBuffer *density_buffer = gpu_buffer_alloc(ctx, grid_size * sizeof(float));

        if (!particle_buffer || !density_buffer) {
            fprintf(stderr, "Failed to allocate GPU buffers: %s\n", gpu_get_error(ctx));
        } else {
            t0 = get_time();
            bool success = gpu_calculate_density(ctx, density_buffer, particle_buffer,
                                                 num_particles, nx, ny, nz, dx,
                                                 origin, origin, origin,
                                                 form_factor_range);
            double gpu_time = get_time() - t0;

            if (success) {
                printf("  GPU time: %.3f ms\n", gpu_time * 1000.0);
                printf("  Speedup: %.2fx\n", cpu_time / gpu_time);

                /* Download results for verification */
                gpu_buffer_download(ctx, density_buffer, density_gpu, grid_size * sizeof(float));

                /* Compare results */
                float max_diff = 0.0f;
                for (size_t i = 0; i < grid_size; i++) {
                    float diff = fabsf(density_cpu[i] - density_gpu[i]);
                    if (diff > max_diff) max_diff = diff;
                }
                printf("  Max CPU-GPU difference: %.6e\n", max_diff);

                if (max_diff < 1e-5f) {
                    printf("  Result: PASS\n");
                } else {
                    printf("  Result: FAIL (values differ)\n");
                }
            } else {
                fprintf(stderr, "GPU density calculation failed: %s\n", gpu_get_error(ctx));
            }

            gpu_buffer_free(particle_buffer);
            gpu_buffer_free(density_buffer);
        }
    }

    /* Cleanup */
    free(particles);
    free(density_cpu);
    free(density_gpu);
    gpu_cleanup(ctx);

    printf("\nTest complete.\n");
    return 0;
}
