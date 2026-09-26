#ifndef LHBUU_UTILS_H
#define LHBUU_UTILS_H

#include <stdio.h>
#include <math.h>

/* Random number generation */
void random_init(unsigned int seed);
float random_uniform(void);  /* Returns [0,1] */
float random_gaussian(float mean, float sigma);

/* Error handling */
void error_exit(const char *fmt, ...);
void warning(const char *fmt, ...);

/* File utilities */
FILE *file_open(const char *filename, const char *mode);

/* Vector operations */
static inline float dot_product(float x1, float y1, float z1, float x2, float y2, float z2) {
    return x1*x2 + y1*y2 + z1*z2;
}

static inline float vector_magnitude(float x, float y, float z) {
    return sqrtf(x*x + y*y + z*z);
}

/* Relativistic energy */
static inline float energy_relativistic(float px, float py, float pz, float mass) {
    float p2 = px*px + py*py + pz*pz;
    return sqrtf(p2 + mass*mass);
}

/* Lorentz transformations */
void lorentz_boost(float *px, float *py, float *pz, float *E,
                   float beta_x, float beta_y, float beta_z);

#endif /* LHBUU_UTILS_H */
