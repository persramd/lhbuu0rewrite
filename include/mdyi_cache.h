/*
 * MDYI Cache: local_pig/g_store System
 *
 * Dynamic per-cell cache for momentum-dependent potential (MDYI mode)
 * Caches which particles contribute to each grid cell and their form factors
 *
 * From original lib/U_alpha_p_md.c lines 30-37, 88-160
 *
 * Key optimization: Form factors depend only on positions (unchanged during
 * force calculation). Cache them once per cell, reuse for all particles
 * evaluated at that cell. ~40x speedup on form factor calculations.
 */

#ifndef MDYI_CACHE_H
#define MDYI_CACHE_H

#include "constants.h"

/*
 * Per-cell cache entry
 *
 * Stores list of particles contributing to a grid cell and their form factors
 */
typedef struct {
    int count;            /* Number of contributing particles */
    int *particle_ids;    /* Array of particle indices [1..count] */
    double *form_factors; /* Array of form factor values [1..count] */
} MDYICacheEntry;

/*
 * MDYI cache structure
 *
 * 3D grid of dynamically-allocated particle lists
 */
typedef struct {
    MDYICacheEntry ***entries;  /* [i][j][k] */
    int dims[3];                 /* Grid dimensions */
    int initialized;             /* First-time setup flag */
} MDYICache;

/*
 * Initialize MDYI cache
 *
 * Args:
 *   nx, ny, nz: grid dimensions
 *
 * Returns:
 *   Allocated MDYICache structure
 */
MDYICache* mdyi_cache_init(int nx, int ny, int nz);

/*
 * Free MDYI cache
 */
void mdyi_cache_free(MDYICache *cache);

/*
 * Clear MDYI cache (free all dynamic lists)
 *
 * Called with signal p_index=0 from force calculation
 * Original: lines 66-76
 *
 * Args:
 *   cache: MDYI cache structure
 */
void mdyi_cache_clear(MDYICache *cache);

/*
 * Check if cell has cached data
 *
 * Args:
 *   cache: MDYI cache structure
 *   i, j, k: grid cell indices
 *
 * Returns:
 *   1 if cached, 0 if not
 */
static inline int mdyi_cache_has_entry(const MDYICache *cache, int i, int j, int k)
{
    return cache->entries[i][j][k].count > 0;
}

/*
 * Get cached particle count
 */
static inline int mdyi_cache_get_count(const MDYICache *cache, int i, int j, int k)
{
    return cache->entries[i][j][k].count;
}

/*
 * Get cached particle ID
 *
 * Args:
 *   n: particle index in list (1-based: 1..count)
 */
static inline int mdyi_cache_get_particle(const MDYICache *cache, int i, int j, int k, int n)
{
    return cache->entries[i][j][k].particle_ids[n];
}

/*
 * Get cached form factor
 *
 * Args:
 *   n: particle index in list (1-based: 1..count)
 */
static inline double mdyi_cache_get_form_factor(const MDYICache *cache, int i, int j, int k, int n)
{
    return cache->entries[i][j][k].form_factors[n];
}

/*
 * Store particle list in cache
 *
 * Called after first evaluation at a grid cell
 * Original: lines 133-143
 *
 * Args:
 *   cache: MDYI cache structure
 *   i, j, k: grid cell indices
 *   count: number of particles
 *   particle_ids: array of particle indices [1..count]
 *   form_factors: array of form factor values [1..count]
 */
void mdyi_cache_store(MDYICache *cache, int i, int j, int k,
                      int count, const int *particle_ids, const double *form_factors);

#endif /* MDYI_CACHE_H */
