#ifndef LHBUU_GPU_COMPUTE_H
#define LHBUU_GPU_COMPUTE_H

#include "types.h"
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* GPU compute backend types */
typedef enum {
    GPU_BACKEND_NONE = 0,    /* CPU fallback */
    GPU_BACKEND_METAL,       /* Apple Metal (macOS) */
    GPU_BACKEND_HIP          /* AMD HIP (Linux/ROCm) - placeholder */
} GPUBackend;

/* GPU device capabilities */
typedef struct {
    GPUBackend backend;
    char device_name[256];
    size_t max_buffer_size;
    size_t max_threads_per_threadgroup;
    size_t total_memory;
    bool unified_memory;     /* True on M-series Macs */
    int compute_units;
} GPUDeviceInfo;

/* Opaque GPU context handle */
typedef struct GPUContext GPUContext;

/* GPU buffer handle for device memory */
typedef struct {
    void *device_ptr;        /* Metal MTLBuffer or HIP device pointer */
    void *host_ptr;          /* CPU-accessible pointer (unified memory) */
    size_t size;
    bool is_unified;         /* Zero-copy access on M-series */
    void *backend_handle;    /* Backend-specific handle */
} GPUBuffer;

/* GPU kernel identifiers */
typedef enum {
    GPU_KERNEL_DENSITY,           /* Density deposition */
    GPU_KERNEL_POTENTIAL,         /* Mean field potential */
    GPU_KERNEL_FORCE,             /* Force calculation */
    GPU_KERNEL_INTEGRATE,         /* Time integration */
    GPU_KERNEL_COULOMB_FFT,       /* Coulomb FFT solver */
    GPU_KERNEL_COUNT
} GPUKernelType;

/* ----- Initialization and Cleanup ----- */

/* Initialize GPU compute context
 * Returns NULL on failure (falls back to CPU) */
GPUContext* gpu_init(void);

/* Get device information */
GPUDeviceInfo gpu_get_device_info(const GPUContext *ctx);

/* Check if GPU is available */
bool gpu_is_available(const GPUContext *ctx);

/* Cleanup GPU context */
void gpu_cleanup(GPUContext *ctx);

/* ----- Memory Management ----- */

/* Allocate GPU buffer (unified memory on M-series) */
GPUBuffer* gpu_buffer_alloc(GPUContext *ctx, size_t size);

/* Allocate and initialize from host data */
GPUBuffer* gpu_buffer_alloc_with_data(GPUContext *ctx, const void *data, size_t size);

/* Free GPU buffer */
void gpu_buffer_free(GPUBuffer *buffer);

/* Copy data to GPU buffer (no-op if unified memory) */
bool gpu_buffer_upload(GPUContext *ctx, GPUBuffer *buffer, const void *data, size_t size);

/* Copy data from GPU buffer (no-op if unified memory) */
bool gpu_buffer_download(GPUContext *ctx, GPUBuffer *buffer, void *data, size_t size);

/* Synchronize buffer (memory barrier on unified memory) */
void gpu_buffer_sync(GPUContext *ctx, GPUBuffer *buffer);

/* ----- Kernel Compilation and Execution ----- */

/* Compile kernel from Metal shader library or HIP source */
bool gpu_kernel_compile(GPUContext *ctx, GPUKernelType kernel);

/* Execute density calculation kernel */
bool gpu_calculate_density(
    GPUContext *ctx,
    GPUBuffer *density_buffer,
    GPUBuffer *particle_buffer,
    size_t num_particles,
    int nx, int ny, int nz,
    float dx,
    float origin_x, float origin_y, float origin_z,
    int form_factor_range
);

/* Execute potential calculation kernel */
bool gpu_calculate_potential(
    GPUContext *ctx,
    GPUBuffer *potential_buffer,
    GPUBuffer *density_buffer,
    int nx, int ny, int nz,
    const PotentialParams *params
);

/* ----- Synchronization ----- */

/* Wait for all GPU operations to complete */
void gpu_synchronize(GPUContext *ctx);

/* ----- Error Handling ----- */

/* Get last error message */
const char* gpu_get_error(const GPUContext *ctx);

/* ----- CPU Fallback Implementations ----- */

/* CPU fallback for density calculation (used when GPU unavailable) */
void cpu_calculate_density(
    float *density,
    const Particle *particles,
    size_t num_particles,
    int nx, int ny, int nz,
    float dx,
    float origin_x, float origin_y, float origin_z,
    int form_factor_range
);

/* CPU fallback for potential calculation */
void cpu_calculate_potential(
    float *potential,
    const float *density,
    int nx, int ny, int nz,
    const PotentialParams *params
);

#ifdef __cplusplus
}
#endif

#endif /* LHBUU_GPU_COMPUTE_H */
