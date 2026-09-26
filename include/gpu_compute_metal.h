#ifndef GPU_COMPUTE_METAL_H
#define GPU_COMPUTE_METAL_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Initialize Metal compute infrastructure */
int gpu_init_metal(void);

/* Cleanup Metal resources */
void gpu_cleanup_metal(void);

/* Grid density calculation using Metal compute shaders */
void grid_calculate_density_metal(Grid *grid, const Particle *particles,
                                   size_t n_particles, int ff_range);

/* Get GPU timing statistics */
void gpu_get_stats_metal(double *total_time, int *call_count);

#ifdef __cplusplus
}
#endif

#endif /* GPU_COMPUTE_METAL_H */
