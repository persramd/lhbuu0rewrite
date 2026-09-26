#include "utils.h"
#include <stdlib.h>
#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

/* PCG Random Number Generator (modern, fast) */
static uint64_t rng_state = 0x4d595df4d0f33173;
static uint64_t rng_inc = 1442695040888963407;

void random_init(unsigned int seed) {
    rng_state = seed;
    rng_inc = (seed << 1) | 1;
}

static uint32_t pcg32(void) {
    uint64_t oldstate = rng_state;
    rng_state = oldstate * 6364136223846793005ULL + rng_inc;
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

float random_uniform(void) {
    return (float)pcg32() / (float)UINT32_MAX;
}

float random_gaussian(float mean, float sigma) {
    /* Box-Muller transform */
    static int have_spare = 0;
    static float spare;

    if (have_spare) {
        have_spare = 0;
        return mean + sigma * spare;
    }

    have_spare = 1;
    float u, v, s;
    do {
        u = random_uniform() * 2.0f - 1.0f;
        v = random_uniform() * 2.0f - 1.0f;
        s = u * u + v * v;
    } while (s >= 1.0f || s == 0.0f);

    s = sqrtf(-2.0f * logf(s) / s);
    spare = v * s;
    return mean + sigma * u * s;
}

void error_exit(const char *fmt, ...) {
    fprintf(stderr, "ERROR: ");
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
    exit(1);
}

void warning(const char *fmt, ...) {
    fprintf(stderr, "WARNING: ");
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
}

FILE *file_open(const char *filename, const char *mode) {
    FILE *fp = fopen(filename, mode);
    if (!fp) {
        error_exit("Cannot open file '%s' with mode '%s'", filename, mode);
    }
    return fp;
}

void lorentz_boost(float *px, float *py, float *pz, float *E,
                   float beta_x, float beta_y, float beta_z) {
    float beta2 = beta_x*beta_x + beta_y*beta_y + beta_z*beta_z;
    if (beta2 == 0.0f) return;

    float gamma = 1.0f / sqrtf(1.0f - beta2);
    float bp = beta_x*(*px) + beta_y*(*py) + beta_z*(*pz);
    float gamma2 = (gamma - 1.0f) / beta2;

    float px_new = *px + gamma2*bp*beta_x + gamma*beta_x*(*E);
    float py_new = *py + gamma2*bp*beta_y + gamma*beta_y*(*E);
    float pz_new = *pz + gamma2*bp*beta_z + gamma*beta_z*(*E);
    float E_new = gamma*(*E + bp);

    *px = px_new;
    *py = py_new;
    *pz = pz_new;
    *E = E_new;
}
