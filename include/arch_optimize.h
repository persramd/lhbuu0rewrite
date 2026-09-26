/*
 * Architecture-Specific Optimizations
 *
 * Detects CPU architecture and provides optimized implementations
 * - M4 Mac: NEON SIMD, 128KB L1, 12MB L2
 * - AMD Ryzen: AVX2 SIMD, 256KB L1, 8MB L2, 16MB L3
 */

#ifndef ARCH_OPTIMIZE_H
#define ARCH_OPTIMIZE_H

#include <stddef.h>
#include <stdint.h>

/* Detect architecture at compile time */
#if defined(__aarch64__) || defined(__arm64__)
    #define ARCH_ARM64 1
    #define CACHE_LINE_SIZE 128  /* M4 cache line */
    #define L1_CACHE_SIZE (128 * 1024)
#elif defined(__x86_64__) || defined(_M_X64)
    #define ARCH_X86_64 1
    #define CACHE_LINE_SIZE 64   /* AMD Ryzen cache line */
    #define L1_CACHE_SIZE (256 * 1024)
#endif

/* Cache-aligned allocation */
#ifdef _WIN32
    #include <malloc.h>
    #define aligned_alloc(align, size) _aligned_malloc(size, align)
    #define aligned_free(ptr) _aligned_free(ptr)
#else
    #include <stdlib.h>
    #define aligned_free(ptr) free(ptr)
#endif

/*
 * Allocate cache-aligned memory
 *
 * Args:
 *   size: number of bytes
 *
 * Returns:
 *   Pointer aligned to CACHE_LINE_SIZE boundary
 */
static inline void* cache_aligned_alloc(size_t size)
{
    return aligned_alloc(CACHE_LINE_SIZE, size);
}

/*
 * SIMD vector dot product (architecture-specific)
 *
 * Computes: result = sum(a[i] * b[i]) for i in [0..n)
 *
 * Args:
 *   a, b: input arrays (must be cache-aligned)
 *   n: number of elements (must be multiple of vector width)
 *
 * Returns:
 *   Dot product
 */
double simd_dot_product(const double *a, const double *b, int n);

/*
 * SIMD vector addition (architecture-specific)
 *
 * Computes: result[i] = a[i] + b[i]
 *
 * Args:
 *   a, b: input arrays (must be cache-aligned)
 *   result: output array (must be cache-aligned)
 *   n: number of elements (must be multiple of vector width)
 */
void simd_vector_add(const double *a, const double *b, double *result, int n);

/*
 * Prefetch data into cache
 *
 * Hints to CPU to load data before it's needed
 *
 * Args:
 *   addr: address to prefetch
 */
static inline void prefetch(const void *addr)
{
#ifdef ARCH_ARM64
    __builtin_prefetch(addr, 0, 3);  /* Read, high temporal locality */
#elif defined(ARCH_X86_64)
    __builtin_prefetch(addr, 0, 3);  /* Read, high temporal locality */
#endif
}

#endif /* ARCH_OPTIMIZE_H */
