#include "gpu_compute.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* CPU-only fallback implementation when GPU is unavailable */

GPUContext* gpu_init(void) {
    fprintf(stderr, "Warning: GPU acceleration not available - using CPU fallback\n");
    return NULL;
}

GPUDeviceInfo gpu_get_device_info(const GPUContext *ctx) {
    (void)ctx;
    GPUDeviceInfo info = {
        .backend = GPU_BACKEND_NONE,
        .device_name = "CPU Fallback",
        .max_buffer_size = 0,
        .max_threads_per_threadgroup = 0,
        .total_memory = 0,
        .unified_memory = false,
        .compute_units = 0
    };
    strncpy(info.device_name, "CPU Fallback", sizeof(info.device_name) - 1);
    return info;
}

bool gpu_is_available(const GPUContext *ctx) {
    (void)ctx;
    return false;
}

void gpu_cleanup(GPUContext *ctx) {
    (void)ctx;
}

GPUBuffer* gpu_buffer_alloc(GPUContext *ctx, size_t size) {
    (void)ctx;
    (void)size;
    return NULL;
}

GPUBuffer* gpu_buffer_alloc_with_data(GPUContext *ctx, const void *data, size_t size) {
    (void)ctx;
    (void)data;
    (void)size;
    return NULL;
}

void gpu_buffer_free(GPUBuffer *buffer) {
    (void)buffer;
}

bool gpu_buffer_upload(GPUContext *ctx, GPUBuffer *buffer, const void *data, size_t size) {
    (void)ctx;
    (void)buffer;
    (void)data;
    (void)size;
    return false;
}

bool gpu_buffer_download(GPUContext *ctx, GPUBuffer *buffer, void *data, size_t size) {
    (void)ctx;
    (void)buffer;
    (void)data;
    (void)size;
    return false;
}

void gpu_buffer_sync(GPUContext *ctx, GPUBuffer *buffer) {
    (void)ctx;
    (void)buffer;
}

bool gpu_kernel_compile(GPUContext *ctx, GPUKernelType kernel) {
    (void)ctx;
    (void)kernel;
    return false;
}

bool gpu_calculate_density(
    GPUContext *ctx,
    GPUBuffer *density_buffer,
    GPUBuffer *particle_buffer,
    size_t num_particles,
    int nx, int ny, int nz,
    float dx,
    float origin_x, float origin_y, float origin_z,
    int form_factor_range
) {
    (void)ctx;
    (void)density_buffer;
    (void)particle_buffer;
    (void)num_particles;
    (void)nx;
    (void)ny;
    (void)nz;
    (void)dx;
    (void)origin_x;
    (void)origin_y;
    (void)origin_z;
    (void)form_factor_range;
    return false;
}

bool gpu_calculate_potential(
    GPUContext *ctx,
    GPUBuffer *potential_buffer,
    GPUBuffer *density_buffer,
    int nx, int ny, int nz,
    const PotentialParams *params
) {
    (void)ctx;
    (void)potential_buffer;
    (void)density_buffer;
    (void)nx;
    (void)ny;
    (void)nz;
    (void)params;
    return false;
}

void gpu_synchronize(GPUContext *ctx) {
    (void)ctx;
}

const char* gpu_get_error(const GPUContext *ctx) {
    (void)ctx;
    return "GPU not available - CPU fallback mode";
}

/* CPU fallback implementations - actual computation code */

void cpu_calculate_density(
    float *density,
    const Particle *particles,
    size_t num_particles,
    int nx, int ny, int nz,
    float dx,
    float origin_x, float origin_y, float origin_z,
    int form_factor_range
) {
    /* Zero out density grid */
    size_t grid_size = nx * ny * nz;
    memset(density, 0, grid_size * sizeof(float));

    /* Simplified form factor: Gaussian deposition */
    const float inv_dx = 1.0f / dx;
    const float sigma_ff = form_factor_range * dx / 2.0f;
    const float norm = 1.0f / (sigma_ff * sqrtf(2.0f * M_PI));

    for (size_t i = 0; i < num_particles; i++) {
        const Particle *p = &particles[i];

        /* Find grid cell */
        int ix = (int)((p->x - origin_x) * inv_dx + 0.5f);
        int iy = (int)((p->y - origin_y) * inv_dx + 0.5f);
        int iz = (int)((p->z - origin_z) * inv_dx + 0.5f);

        /* Deposit density in neighborhood */
        for (int dx_cell = -form_factor_range; dx_cell <= form_factor_range; dx_cell++) {
            int gx = ix + dx_cell;
            if (gx < 0 || gx >= nx) continue;

            for (int dy_cell = -form_factor_range; dy_cell <= form_factor_range; dy_cell++) {
                int gy = iy + dy_cell;
                if (gy < 0 || gy >= ny) continue;

                for (int dz_cell = -form_factor_range; dz_cell <= form_factor_range; dz_cell++) {
                    int gz = iz + dz_cell;
                    if (gz < 0 || gz >= nz) continue;

                    /* Distance from particle to cell center */
                    float cell_x = origin_x + gx * dx;
                    float cell_y = origin_y + gy * dx;
                    float cell_z = origin_z + gz * dx;

                    float dist2 = (p->x - cell_x) * (p->x - cell_x) +
                                  (p->y - cell_y) * (p->y - cell_y) +
                                  (p->z - cell_z) * (p->z - cell_z);

                    /* Gaussian weight */
                    float weight = norm * expf(-0.5f * dist2 / (sigma_ff * sigma_ff));

                    int idx = gx * ny * nz + gy * nz + gz;
                    density[idx] += weight;
                }
            }
        }
    }

    /* Normalize by cell volume */
    float cell_volume = dx * dx * dx;
    for (size_t i = 0; i < grid_size; i++) {
        density[i] /= cell_volume;
    }
}

void cpu_calculate_potential(
    float *potential,
    const float *density,
    int nx, int ny, int nz,
    const PotentialParams *params
) {
    size_t grid_size = nx * ny * nz;

    /* Skyrme potential: U = A*ρ + B*ρ^σ */
    for (size_t i = 0; i < grid_size; i++) {
        float rho = density[i];
        potential[i] = params->A * rho +
                       params->B * powf(rho, params->sigma);
    }
}
