#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include "gpu_compute.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Metal GPU context */
struct GPUContext {
    id<MTLDevice> device;
    id<MTLCommandQueue> command_queue;
    id<MTLLibrary> library;
    id<MTLComputePipelineState> pipelines[GPU_KERNEL_COUNT];
    char error_msg[512];
    GPUDeviceInfo info;
    bool initialized;
};

/* GPU buffer wrapper */
struct GPUBufferImpl {
    id<MTLBuffer> metal_buffer;
    void *host_ptr;
    size_t size;
    bool is_unified;
};

/* ----- Initialization ----- */

GPUContext* gpu_init(void) {
    GPUContext *ctx = calloc(1, sizeof(GPUContext));
    if (!ctx) {
        return NULL;
    }

    /* Get default Metal device (M4 GPU) */
    ctx->device = MTLCreateSystemDefaultDevice();
    if (!ctx->device) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Metal device not found - GPU unavailable");
        free(ctx);
        return NULL;
    }

    /* Create command queue for GPU work submission */
    ctx->command_queue = [ctx->device newCommandQueue];
    if (!ctx->command_queue) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Failed to create Metal command queue");
        free(ctx);
        return NULL;
    }

    /* Query device capabilities */
    ctx->info.backend = GPU_BACKEND_METAL;
    strncpy(ctx->info.device_name,
            [[ctx->device name] UTF8String],
            sizeof(ctx->info.device_name) - 1);

    /* M4 has unified memory architecture */
    ctx->info.unified_memory = true;
    ctx->info.max_buffer_size = [ctx->device maxBufferLength];
    ctx->info.max_threads_per_threadgroup = [ctx->device maxThreadsPerThreadgroup].width;

    /* M4 Mac: 10-core GPU on base model, 16GB-32GB unified memory */
    if ([ctx->device supportsFamily:MTLGPUFamilyApple9]) {
        /* Apple9 = M4 generation */
        ctx->info.compute_units = 10;  /* 10-core GPU on M4 base */
    } else if ([ctx->device supportsFamily:MTLGPUFamilyApple8]) {
        /* Apple8 = M3 generation */
        ctx->info.compute_units = 8;
    } else {
        /* Older Apple Silicon */
        ctx->info.compute_units = 7;
    }

    /* Get recommended working set size (unified memory pool) */
    ctx->info.total_memory = [ctx->device recommendedMaxWorkingSetSize];

    ctx->initialized = true;
    snprintf(ctx->error_msg, sizeof(ctx->error_msg), "No error");

    printf("Metal GPU initialized: %s\n", ctx->info.device_name);
    printf("  Compute units: %d cores\n", ctx->info.compute_units);
    printf("  Unified memory: %zu MB\n", ctx->info.total_memory / (1024*1024));
    printf("  Max threads/threadgroup: %zu\n", ctx->info.max_threads_per_threadgroup);

    return ctx;
}

GPUDeviceInfo gpu_get_device_info(const GPUContext *ctx) {
    if (ctx) {
        return ctx->info;
    }
    GPUDeviceInfo empty = {0};
    return empty;
}

bool gpu_is_available(const GPUContext *ctx) {
    return ctx && ctx->initialized && ctx->device != nil;
}

void gpu_cleanup(GPUContext *ctx) {
    if (!ctx) return;

    /* Release compute pipelines */
    for (int i = 0; i < GPU_KERNEL_COUNT; i++) {
        if (ctx->pipelines[i]) {
            ctx->pipelines[i] = nil;
        }
    }

    /* Release Metal objects (ARC handles this) */
    ctx->library = nil;
    ctx->command_queue = nil;
    ctx->device = nil;

    free(ctx);
}

/* ----- Memory Management ----- */

GPUBuffer* gpu_buffer_alloc(GPUContext *ctx, size_t size) {
    if (!ctx || !ctx->device) {
        return NULL;
    }

    GPUBuffer *buffer = calloc(1, sizeof(GPUBuffer));
    if (!buffer) {
        return NULL;
    }

    /* Allocate unified memory buffer (zero-copy on M4)
     * MTLResourceStorageModeShared: CPU and GPU share same memory */
    id<MTLBuffer> mtl_buffer = [ctx->device newBufferWithLength:size
                                                        options:MTLResourceStorageModeShared];

    if (!mtl_buffer) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Failed to allocate Metal buffer of size %zu", size);
        free(buffer);
        return NULL;
    }

    buffer->device_ptr = (__bridge void*)mtl_buffer;
    buffer->host_ptr = [mtl_buffer contents];  /* Direct CPU access */
    buffer->size = size;
    buffer->is_unified = true;
    buffer->backend_handle = (__bridge_retained void*)mtl_buffer;

    return buffer;
}

GPUBuffer* gpu_buffer_alloc_with_data(GPUContext *ctx, const void *data, size_t size) {
    GPUBuffer *buffer = gpu_buffer_alloc(ctx, size);
    if (!buffer) {
        return NULL;
    }

    /* Copy data directly to unified memory */
    memcpy(buffer->host_ptr, data, size);

    return buffer;
}

void gpu_buffer_free(GPUBuffer *buffer) {
    if (!buffer) return;

    if (buffer->backend_handle) {
        /* Release Metal buffer */
        CFRelease(buffer->backend_handle);
    }

    free(buffer);
}

bool gpu_buffer_upload(GPUContext *ctx, GPUBuffer *buffer, const void *data, size_t size) {
    if (!ctx || !buffer || !data) {
        return false;
    }

    if (size > buffer->size) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Upload size %zu exceeds buffer size %zu", size, buffer->size);
        return false;
    }

    /* On unified memory, just memcpy (zero-copy architecture) */
    memcpy(buffer->host_ptr, data, size);

    return true;
}

bool gpu_buffer_download(GPUContext *ctx, GPUBuffer *buffer, void *data, size_t size) {
    if (!ctx || !buffer || !data) {
        return false;
    }

    if (size > buffer->size) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Download size %zu exceeds buffer size %zu", size, buffer->size);
        return false;
    }

    /* On unified memory, just memcpy */
    memcpy(data, buffer->host_ptr, size);

    return true;
}

void gpu_buffer_sync(GPUContext *ctx, GPUBuffer *buffer) {
    if (!ctx || !buffer) return;

    /* On unified memory, Metal handles coherency automatically
     * No explicit synchronization needed for StorageModeShared */
    (void)ctx;
    (void)buffer;
}

/* ----- Kernel Compilation ----- */

bool gpu_kernel_compile(GPUContext *ctx, GPUKernelType kernel) {
    if (!ctx || !ctx->device) {
        return false;
    }

    /* Load default Metal library (compiled .metallib from .metal shaders)
     * Metal shaders will be in separate .metal files compiled at build time */
    NSError *error = nil;

    if (!ctx->library) {
        /* Try loading from default library (built into app bundle) */
        ctx->library = [ctx->device newDefaultLibrary];

        if (!ctx->library) {
            snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                     "Failed to load Metal shader library - shaders not compiled");
            return false;
        }
    }

    /* Get kernel function by name */
    const char *kernel_names[] = {
        "calculate_density",
        "calculate_potential",
        "calculate_force",
        "integrate_particles",
        "coulomb_fft"
    };

    if (kernel < 0 || kernel >= GPU_KERNEL_COUNT) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Invalid kernel type %d", kernel);
        return false;
    }

    id<MTLFunction> function = [ctx->library
        newFunctionWithName:[NSString stringWithUTF8String:kernel_names[kernel]]];

    if (!function) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Kernel function '%s' not found in Metal library", kernel_names[kernel]);
        return false;
    }

    /* Create compute pipeline state */
    ctx->pipelines[kernel] = [ctx->device newComputePipelineStateWithFunction:function
                                                                         error:&error];

    if (!ctx->pipelines[kernel]) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Failed to compile kernel '%s': %s",
                 kernel_names[kernel],
                 [[error localizedDescription] UTF8String]);
        return false;
    }

    return true;
}

/* ----- Kernel Execution ----- */

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
    if (!ctx || !ctx->device || !density_buffer || !particle_buffer) {
        return false;
    }

    /* Ensure kernel is compiled */
    if (!ctx->pipelines[GPU_KERNEL_DENSITY]) {
        if (!gpu_kernel_compile(ctx, GPU_KERNEL_DENSITY)) {
            return false;
        }
    }

    /* Create command buffer */
    id<MTLCommandBuffer> command_buffer = [ctx->command_queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];

    /* Set pipeline and buffers */
    [encoder setComputePipelineState:ctx->pipelines[GPU_KERNEL_DENSITY]];

    id<MTLBuffer> density_mtl = (__bridge id<MTLBuffer>)density_buffer->device_ptr;
    id<MTLBuffer> particle_mtl = (__bridge id<MTLBuffer>)particle_buffer->device_ptr;

    [encoder setBuffer:density_mtl offset:0 atIndex:0];
    [encoder setBuffer:particle_mtl offset:0 atIndex:1];

    /* Set kernel parameters (passed as uniform buffer or individual values) */
    struct {
        uint32_t num_particles;
        int32_t nx, ny, nz;
        float dx;
        float origin_x, origin_y, origin_z;
        int32_t form_factor_range;
    } params = {
        .num_particles = (uint32_t)num_particles,
        .nx = nx, .ny = ny, .nz = nz,
        .dx = dx,
        .origin_x = origin_x, .origin_y = origin_y, .origin_z = origin_z,
        .form_factor_range = form_factor_range
    };

    [encoder setBytes:&params length:sizeof(params) atIndex:2];

    /* Dispatch kernel: 1D grid over particles */
    MTLSize grid_size = MTLSizeMake(num_particles, 1, 1);
    NSUInteger threads_per_group = ctx->pipelines[GPU_KERNEL_DENSITY].maxTotalThreadsPerThreadgroup;
    MTLSize threadgroup_size = MTLSizeMake(MIN(threads_per_group, num_particles), 1, 1);

    [encoder dispatchThreads:grid_size threadsPerThreadgroup:threadgroup_size];
    [encoder endEncoding];

    /* Submit and wait */
    [command_buffer commit];
    [command_buffer waitUntilCompleted];

    if (command_buffer.error) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Density kernel execution failed: %s",
                 [[command_buffer.error localizedDescription] UTF8String]);
        return false;
    }

    return true;
}

bool gpu_calculate_potential(
    GPUContext *ctx,
    GPUBuffer *potential_buffer,
    GPUBuffer *density_buffer,
    int nx, int ny, int nz,
    const PotentialParams *params
) {
    if (!ctx || !ctx->device || !potential_buffer || !density_buffer || !params) {
        return false;
    }

    /* Ensure kernel is compiled */
    if (!ctx->pipelines[GPU_KERNEL_POTENTIAL]) {
        if (!gpu_kernel_compile(ctx, GPU_KERNEL_POTENTIAL)) {
            return false;
        }
    }

    /* Create command buffer */
    id<MTLCommandBuffer> command_buffer = [ctx->command_queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];

    /* Set pipeline and buffers */
    [encoder setComputePipelineState:ctx->pipelines[GPU_KERNEL_POTENTIAL]];

    id<MTLBuffer> potential_mtl = (__bridge id<MTLBuffer>)potential_buffer->device_ptr;
    id<MTLBuffer> density_mtl = (__bridge id<MTLBuffer>)density_buffer->device_ptr;

    [encoder setBuffer:potential_mtl offset:0 atIndex:0];
    [encoder setBuffer:density_mtl offset:0 atIndex:1];

    /* Set grid dimensions and potential parameters */
    struct {
        int32_t nx, ny, nz;
        float A, B, sigma;
        float C, lambda;
        float D;
        float A_surf, lambda_surf;
        float rho0;
    } kernel_params = {
        .nx = nx, .ny = ny, .nz = nz,
        .A = params->A, .B = params->B, .sigma = params->sigma,
        .C = params->C, .lambda = params->lambda,
        .D = params->D,
        .A_surf = params->A_surf, .lambda_surf = params->lambda_surf,
        .rho0 = params->rho0
    };

    [encoder setBytes:&kernel_params length:sizeof(kernel_params) atIndex:2];

    /* Dispatch kernel: 3D grid over cells */
    size_t total_cells = nx * ny * nz;
    MTLSize grid_size = MTLSizeMake(total_cells, 1, 1);
    NSUInteger threads_per_group = ctx->pipelines[GPU_KERNEL_POTENTIAL].maxTotalThreadsPerThreadgroup;
    MTLSize threadgroup_size = MTLSizeMake(MIN(threads_per_group, total_cells), 1, 1);

    [encoder dispatchThreads:grid_size threadsPerThreadgroup:threadgroup_size];
    [encoder endEncoding];

    /* Submit and wait */
    [command_buffer commit];
    [command_buffer waitUntilCompleted];

    if (command_buffer.error) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg),
                 "Potential kernel execution failed: %s",
                 [[command_buffer.error localizedDescription] UTF8String]);
        return false;
    }

    return true;
}

/* ----- Synchronization ----- */

void gpu_synchronize(GPUContext *ctx) {
    if (!ctx || !ctx->command_queue) return;

    /* Create and immediately wait on empty command buffer */
    id<MTLCommandBuffer> sync_buffer = [ctx->command_queue commandBuffer];
    [sync_buffer commit];
    [sync_buffer waitUntilCompleted];
}

/* ----- Error Handling ----- */

const char* gpu_get_error(const GPUContext *ctx) {
    if (!ctx) {
        return "GPU context is NULL";
    }
    return ctx->error_msg;
}

/* ----- CPU Fallback Implementations ----- */

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

    /* Simplified form factor: Gaussian deposition
     * Real implementation would use Woods-Saxon or similar */
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
