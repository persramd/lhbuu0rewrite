/*
 * MDYI Cache Implementation
 *
 * Dynamic per-cell caching for momentum-dependent potential
 * Direct port from lib/U_alpha_p_md.c
 */

#include "mdyi_cache.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 * Initialize MDYI cache
 */
MDYICache* mdyi_cache_init(int nx, int ny, int nz)
{
    MDYICache *cache = malloc(sizeof(MDYICache));

    cache->dims[0] = nx;
    cache->dims[1] = ny;
    cache->dims[2] = nz;
    cache->initialized = 1;

    /* Allocate 3D grid of cache entries (1-based indexing) */
    cache->entries = malloc((nx + 1) * sizeof(MDYICacheEntry**));

    for (int i = 0; i <= nx; i++) {
        cache->entries[i] = malloc((ny + 1) * sizeof(MDYICacheEntry*));

        for (int j = 0; j <= ny; j++) {
            cache->entries[i][j] = calloc(nz + 1, sizeof(MDYICacheEntry));

            /* Initialize all entries as empty */
            for (int k = 0; k <= nz; k++) {
                cache->entries[i][j][k].count = 0;
                cache->entries[i][j][k].particle_ids = NULL;
                cache->entries[i][j][k].form_factors = NULL;
            }
        }
    }

    return cache;
}

/*
 * Free MDYI cache
 */
void mdyi_cache_free(MDYICache *cache)
{
    if (!cache) return;

    /* Free all dynamic arrays first */
    for (int i = 0; i <= cache->dims[0]; i++) {
        for (int j = 0; j <= cache->dims[1]; j++) {
            for (int k = 0; k <= cache->dims[2]; k++) {
                if (cache->entries[i][j][k].count > 0) {
                    free(cache->entries[i][j][k].particle_ids);
                    free(cache->entries[i][j][k].form_factors);
                }
            }
            free(cache->entries[i][j]);
        }
        free(cache->entries[i]);
    }
    free(cache->entries);
    free(cache);
}

/*
 * Clear MDYI cache
 *
 * Direct port from lines 66-76
 */
void mdyi_cache_clear(MDYICache *cache)
{
    for (int i = 1; i <= cache->dims[0]; i++) {
        for (int j = 1; j <= cache->dims[1]; j++) {
            for (int k = 1; k <= cache->dims[2]; k++) {
                if (cache->entries[i][j][k].count > 0) {
                    /* Free dynamic arrays */
                    free(cache->entries[i][j][k].particle_ids);
                    free(cache->entries[i][j][k].form_factors);

                    /* Reset entry */
                    cache->entries[i][j][k].count = 0;
                    cache->entries[i][j][k].particle_ids = NULL;
                    cache->entries[i][j][k].form_factors = NULL;
                }
            }
        }
    }
}

/*
 * Store particle list in cache
 *
 * Direct port from lines 133-143
 */
void mdyi_cache_store(MDYICache *cache, int i, int j, int k,
                      int count, const int *particle_ids, const double *form_factors)
{
    if (count <= 0) return;

    MDYICacheEntry *entry = &cache->entries[i][j][k];

    /* Allocate arrays for particle list (1-based indexing, so +1) */
    entry->particle_ids = malloc((count + 1) * sizeof(int));
    entry->form_factors = malloc((count + 1) * sizeof(double));

    if (!entry->particle_ids || !entry->form_factors) {
        fprintf(stderr, "ERROR: Failed to allocate MDYI cache at [%d][%d][%d]\n", i, j, k);
        exit(1);
    }

    /* Copy data (lines 140-142) */
    entry->count = count;
    for (int n = 1; n <= count; n++) {
        entry->particle_ids[n] = particle_ids[n];
        entry->form_factors[n] = form_factors[n];
    }
}
