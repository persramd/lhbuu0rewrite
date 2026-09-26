/*
 * LHBUU Mathematical Macros
 *
 * Philosophy: Change once, applies everywhere
 * - Single point of maintenance
 * - Compiler optimization friendly
 * - Architecture-specific tuning possible
 */

#ifndef MACROS_H
#define MACROS_H

#include <math.h>

/* Basic powers */
#define square(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

/* Vector magnitudes */
#define mag3(x,y,z) sqrt(square(x) + square(y) + square(z))
#define mag4(px,py,pz,m) sqrt(square(px) + square(py) + square(pz) + square(m))

/* Roots */
#define cbrt(x) pow((x), 1.0/3.0)
#define sixth_root(x) pow((x), 1.0/6.0)

/* Geometry */
#define sphere_volume(r) ((4.0*M_PI/3.0) * cube(r))

/* Grid operations */
#define where_am_I(x,g_cntr,inv_lat_spac) ((int)(((x) + (g_cntr)) * (inv_lat_spac)))

/* Sign function */
#define sgn(x) (((x) > 0) ? 1 : -1)

#endif /* MACROS_H */
