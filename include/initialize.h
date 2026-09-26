#ifndef LHBUU_INITIALIZE_H
#define LHBUU_INITIALIZE_H

#include "types.h"

/* Initialize simulation state */
int sim_init(SimState *state, SimConfig *config);

/* Free simulation state */
void sim_free(SimState *state);

/* Generate initial nuclear configuration */
int initialize_nuclei(SimState *state);

/* Load nuclear density profile */
int initialize_load_profile(Nucleus *nucleus, const char *filename);

/* Generate Woods-Saxon profile */
void initialize_woods_saxon(Nucleus *nucleus);

#endif /* LHBUU_INITIALIZE_H */
