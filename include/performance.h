#ifndef LHBUU_PERFORMANCE_H
#define LHBUU_PERFORMANCE_H

#include "types.h"
#include <stddef.h>

/* Performance modes based on working set size */
typedef enum {
    PERF_CACHE_FIT,      /* Everything fits in L2 - simple algorithms */
    PERF_CACHE_AWARE,    /* Fits in L3 - periodic sorting, neighbor lists */
    PERF_MEMORY_BOUND    /* Large system - aggressive optimization */
} PerformanceMode;

/* Memory and cache statistics */
typedef struct {
    size_t particle_memory;
    size_t grid_memory;
    size_t neighbor_memory;
    size_t total_memory;
    size_t l2_cache_size;
    size_t l3_cache_size;
    PerformanceMode mode;
} MemoryProfile;

/* Analyze memory requirements and choose strategy */
PerformanceMode perf_analyze_memory(const SimConfig *config, MemoryProfile *profile);

/* Print memory and performance info */
void perf_print_profile(const MemoryProfile *profile);

/* Get cache sizes (platform-specific) */
size_t perf_get_l2_cache_size(void);
size_t perf_get_l3_cache_size(void);

/* Get number of CPU cores */
int perf_get_num_cores(void);

#endif /* LHBUU_PERFORMANCE_H */
