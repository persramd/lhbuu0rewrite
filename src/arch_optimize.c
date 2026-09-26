/*
 * Architecture-Specific SIMD Implementations
 *
 * NEON for ARM64 (M4 Mac)
 * AVX2 for x86_64 (AMD Ryzen)
 */

#include "arch_optimize.h"

#ifdef ARCH_ARM64
    #include <arm_neon.h>
#elif defined(ARCH_X86_64)
    #include <immintrin.h>
#endif

/*
 * SIMD dot product
 */
double simd_dot_product(const double *a, const double *b, int n)
{
#ifdef ARCH_ARM64
    /* ARM NEON: 2 doubles per vector */
    float64x2_t sum_vec = vdupq_n_f64(0.0);

    for (int i = 0; i < n; i += 2) {
        float64x2_t a_vec = vld1q_f64(&a[i]);
        float64x2_t b_vec = vld1q_f64(&b[i]);
        sum_vec = vfmaq_f64(sum_vec, a_vec, b_vec);  /* FMA: sum += a * b */
    }

    /* Horizontal reduction */
    return vgetq_lane_f64(sum_vec, 0) + vgetq_lane_f64(sum_vec, 1);

#elif defined(ARCH_X86_64)
    /* x86 AVX2: 4 doubles per vector */
    __m256d sum_vec = _mm256_setzero_pd();

    for (int i = 0; i < n; i += 4) {
        __m256d a_vec = _mm256_load_pd(&a[i]);
        __m256d b_vec = _mm256_load_pd(&b[i]);
        sum_vec = _mm256_fmadd_pd(a_vec, b_vec, sum_vec);  /* FMA: sum += a * b */
    }

    /* Horizontal reduction */
    __m128d sum_high = _mm256_extractf128_pd(sum_vec, 1);
    __m128d sum_low = _mm256_castpd256_pd128(sum_vec);
    __m128d sum_128 = _mm_add_pd(sum_low, sum_high);
    __m128d sum_64 = _mm_hadd_pd(sum_128, sum_128);

    return _mm_cvtsd_f64(sum_64);

#else
    /* Fallback: scalar code */
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
#endif
}

/*
 * SIMD vector addition
 */
void simd_vector_add(const double *a, const double *b, double *result, int n)
{
#ifdef ARCH_ARM64
    /* ARM NEON: 2 doubles per vector */
    for (int i = 0; i < n; i += 2) {
        float64x2_t a_vec = vld1q_f64(&a[i]);
        float64x2_t b_vec = vld1q_f64(&b[i]);
        float64x2_t sum_vec = vaddq_f64(a_vec, b_vec);
        vst1q_f64(&result[i], sum_vec);
    }

#elif defined(ARCH_X86_64)
    /* x86 AVX2: 4 doubles per vector */
    for (int i = 0; i < n; i += 4) {
        __m256d a_vec = _mm256_load_pd(&a[i]);
        __m256d b_vec = _mm256_load_pd(&b[i]);
        __m256d sum_vec = _mm256_add_pd(a_vec, b_vec);
        _mm256_store_pd(&result[i], sum_vec);
    }

#else
    /* Fallback: scalar code */
    for (int i = 0; i < n; i++) {
        result[i] = a[i] + b[i];
    }
#endif
}
