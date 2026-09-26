#ifndef LHBUU_DIAGNOSTICS_H
#define LHBUU_DIAGNOSTICS_H

#include "types.h"
#include <stdio.h>

/* Energy calculations */
float diagnostics_kinetic_energy(const SimState *state);
float diagnostics_potential_energy(const SimState *state);
float diagnostics_total_energy(const SimState *state);

/* Momentum conservation */
void diagnostics_momentum(const SimState *state, float *px, float *py, float *pz);

/* Angular momentum */
void diagnostics_angular_momentum(const SimState *state, float *lx, float *ly, float *lz);

/* Write timestep diagnostics */
void diagnostics_write_timestep(FILE *fp, const SimState *state);

/* Write phase space output */
int diagnostics_write_phase_space(const char *filename, const SimState *state);

/* Print summary */
void diagnostics_print_summary(const SimState *state);

#endif /* LHBUU_DIAGNOSTICS_H */
