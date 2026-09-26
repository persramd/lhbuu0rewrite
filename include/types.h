#ifndef LHBUU_TYPES_H
#define LHBUU_TYPES_H

#include <stddef.h>
#include <stdint.h>

/* Physical constants */
#define PI 3.141592653589793
#define HBARC 0.197327  /* GeV·fm */
#define NUCLEON_MASS 0.938919  /* GeV */

/* Particle state */
typedef struct {
    float x, y, z;          /* Position (fm) */
    float px, py, pz;       /* Momentum (GeV) */
    float dxdt, dydt, dzdt; /* Position derivatives */
    float dpxdt, dpydt, dpzdt; /* Momentum derivatives */
    int charge;             /* 0=neutron, 1=proton */
    int ensemble_id;        /* Ensemble number */
    int nucleon_id;         /* ID within nucleus */
    int is_projectile;      /* 1=projectile, 0=target */
    int num_collisions;     /* Collision counter */
    int grid_cell;          /* Current grid cell index */
} Particle;

/* Grid system */
typedef struct {
    int nx, ny, nz;         /* Grid dimensions */
    float dx;               /* Lattice spacing (fm) */
    float origin_x, origin_y, origin_z; /* Grid origin */

    float *baryon_density;  /* ρ(x,y,z) */
    float *charge_density;  /* ρ_charge(x,y,z) */
    float *potential;       /* U(x,y,z) mean field */
    float *coulomb_potential; /* V_coulomb(x,y,z) */

    int **cell_particles;   /* Particle indices per cell */
    int *cell_counts;       /* Particle count per cell */
    int max_per_cell;       /* Max particles per cell */
} Grid;

/* Nuclear configuration */
typedef struct {
    int Z;                  /* Atomic number */
    int A;                  /* Mass number */
    float radius;           /* Nuclear radius (fm) */
    float *density_profile; /* r-dependent density */
} Nucleus;

/* Potential parameters */
typedef struct {
    float A, B, sigma;      /* Skyrme: U = A*ρ + B*ρ^σ */
    float C, lambda;        /* Momentum-dependent */
    float D;                /* Isospin */
    float A_surf, lambda_surf; /* Surface Yukawa */
    float rho0;             /* Saturation density (fm^-3) */
    int enable_coulomb;     /* Boolean */
    int enable_momentum_dep; /* Boolean */
} PotentialParams;

/* Collision parameters */
typedef struct {
    float sigma_nn_max;     /* Max NN cross section (mb) */
    float sigma_rho_reduction; /* Density reduction factor */
    int enable_pauli_blocking; /* Boolean */
    int pauli_ensemble_count; /* Ensembles for Pauli check */
} CollisionParams;

/* Simulation configuration */
typedef struct {
    /* Nuclear system */
    Nucleus projectile;
    Nucleus target;

    /* Kinematics */
    float beam_energy;      /* MeV per nucleon */
    float impact_parameter; /* Normalized (0-1) */

    /* Time evolution */
    float dt;               /* Time step (fm/c) */
    float t_max;            /* Max simulation time (fm/c) */
    int diagnostic_interval; /* Steps between diagnostics */

    /* Ensembles */
    int num_ensembles;      /* Event count */

    /* Grid */
    Grid grid;
    int form_factor_range;  /* Cells */
    int sort_interval;      /* Grid re-sort interval */

    /* Physics */
    PotentialParams potential;
    CollisionParams collision;

    /* Output */
    char output_prefix[256];
    char output_path[256];
    int full_output;        /* Write final phase space */
    int animation_interval; /* Animation frames */
} SimConfig;

/* Simulation state */
typedef struct {
    SimConfig *config;
    Particle *particles;
    size_t num_particles;
    size_t max_particles;
    int current_timestep;
    float current_time;

    /* Statistics */
    int total_collisions;
    int pauli_blocked;
} SimState;

#endif /* LHBUU_TYPES_H */
