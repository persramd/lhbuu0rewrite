#ifndef LHBUU_CONFIG_H
#define LHBUU_CONFIG_H

#include <stdio.h>
#include "types.h"

/* Load configuration from file */
int config_load(const char *filename, SimConfig *config);

/* Set default configuration */
void config_set_defaults(SimConfig *config);

/* Validate configuration */
int config_validate(const SimConfig *config);

/* Print configuration to file */
void config_print(FILE *fp, const SimConfig *config);

#endif /* LHBUU_CONFIG_H */
