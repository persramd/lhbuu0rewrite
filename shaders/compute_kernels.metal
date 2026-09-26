#include <metal_stdlib>
using namespace metal;

/*
 * LHBUU Metal Compute Kernels
 * M4 Mac - Unified Memory Architecture
 */

/* Particle structure matching types.h */
struct Particle {
    float x, y, z;
    float px, py, pz;
    float dxdt, dydt, dzdt;
    float dpxdt, dpydt, dpzdt;
    int charge;
    int ensemble_id;
    int nucleon_id;
    int is_projectile;
    int num_collisions;
    int grid_cell;
};

/* Kernel parameters passed from CPU */
struct DensityParams {
    uint num_particles;
    int nx, ny, nz;
    float dx;
    float origin_x, origin_y, origin_z;
    int form_factor_range;
};

struct PotentialParams {
    int nx, ny, nz;
    float A, B, sigma;
    float C, lambda;
    float D;
    float A_surf, lambda_surf;
    float rho0;
};

/* Form factor (Monaghan cubic spline kernel) - matches grid.c implementation */
static float form_factor_1d(float r) {
    float q = fabs(r);
    if (q >= 2.0f) return 0.0f;
    if (q >= 1.0f) {
        float tmp = 2.0f - q;
        return 0.25f * tmp * tmp * tmp;
    }
    float q2 = q * q;
    return 1.0f - 1.5f * q2 + 0.75f * q2 * q;
}

/* Grid indexing - matches grid.h */
static inline int grid_index(int ix, int iy, int iz, int ny, int nz) {
    return ix * ny * nz + iy * nz + iz;
}

/* Convert world coordinates to grid cell indices - matches grid.h */
static inline void grid_coords(float x, float y, float z,
                                float origin_x, float origin_y, float origin_z,
                                float dx,
                                thread int* ix, thread int* iy, thread int* iz) {
    *ix = (int)((x - origin_x) / dx + 0.5f);
    *iy = (int)((y - origin_y) / dx + 0.5f);
    *iz = (int)((z - origin_z) / dx + 0.5f);
}

/*
 * Density calculation kernel
 *
 * Algorithm (matching src/grid.c::grid_calculate_density):
 * 1. One thread per particle
 * 2. Find particle's home grid cell
 * 3. Loop over ±ff_range neighborhood (typically ±2 cells)
 * 4. Evaluate separable 3D form factor = f_x * f_y * f_z
 * 5. Atomically accumulate to baryon_density
 * 6. If charged, atomically accumulate to charge_density
 *
 * Threading: 1D dispatch over particles
 * Memory: Atomic adds to grid (scattered writes)
 */
kernel void calculate_density(
    device atomic_float* baryon_density [[buffer(0)]],
    device atomic_float* charge_density [[buffer(1)]],
    device const Particle* particles [[buffer(2)]],
    constant DensityParams& params [[buffer(3)]],
    uint particle_idx [[thread_position_in_grid]])
{
    /* Bounds check */
    if (particle_idx >= params.num_particles) return;

    const device Particle& p = particles[particle_idx];

    /* Find particle's home grid cell */
    int ix, iy, iz;
    grid_coords(p.x, p.y, p.z,
                params.origin_x, params.origin_y, params.origin_z,
                params.dx, &ix, &iy, &iz);

    /* Distribute density over neighboring cells within form factor support */
    for (int dx = -params.form_factor_range; dx <= params.form_factor_range; dx++) {
        float ff_x = form_factor_1d((float)dx);
        if (ff_x == 0.0f) continue;

        for (int dy = -params.form_factor_range; dy <= params.form_factor_range; dy++) {
            float ff_y = form_factor_1d((float)dy);
            if (ff_y == 0.0f) continue;

            for (int dz = -params.form_factor_range; dz <= params.form_factor_range; dz++) {
                float ff_z = form_factor_1d((float)dz);
                if (ff_z == 0.0f) continue;

                int cx = ix + dx;
                int cy = iy + dy;
                int cz = iz + dz;

                /* Boundary check (no periodic BC) */
                if (cx < 0 || cx >= params.nx) continue;
                if (cy < 0 || cy >= params.ny) continue;
                if (cz < 0 || cz >= params.nz) continue;

                /* Separable 3D form factor */
                float ff = ff_x * ff_y * ff_z;

                int idx = grid_index(cx, cy, cz, params.ny, params.nz);

                /* Atomic accumulation (thread-safe) */
                atomic_fetch_add_explicit(&baryon_density[idx], ff,
                                         memory_order_relaxed);
                if (p.charge) {
                    atomic_fetch_add_explicit(&charge_density[idx], ff,
                                             memory_order_relaxed);
                }
            }
        }
    }
}

/*
 * Potential calculation kernel
 *
 * Skyrme mean field potential: U = A*ρ + B*ρ^σ
 *
 * Threading: 1D dispatch over grid cells
 * Memory: Read density, write potential (1:1 mapping)
 */
kernel void calculate_potential(
    device float* potential [[buffer(0)]],
    device const float* density [[buffer(1)]],
    constant PotentialParams& params [[buffer(2)]],
    uint cell_idx [[thread_position_in_grid]])
{
    uint n_cells = params.nx * params.ny * params.nz;
    if (cell_idx >= n_cells) return;

    float rho = density[cell_idx];

    /* Skyrme potential: U = A*ρ + B*ρ^σ */
    potential[cell_idx] = params.A * rho + params.B * pow(rho, params.sigma);
}

/*
 * Grid zeroing kernel
 *
 * Zeros density grids before accumulation
 * Separate kernel for better performance than memset from CPU
 */
kernel void zero_density(
    device atomic_float* baryon_density [[buffer(0)]],
    device atomic_float* charge_density [[buffer(1)]],
    constant DensityParams& params [[buffer(2)]],
    uint cell_idx [[thread_position_in_grid]])
{
    uint n_cells = params.nx * params.ny * params.nz;
    if (cell_idx >= n_cells) return;

    atomic_store_explicit(&baryon_density[cell_idx], 0.0f,
                         memory_order_relaxed);
    atomic_store_explicit(&charge_density[cell_idx], 0.0f,
                         memory_order_relaxed);
}

/*
 * Density normalization kernel
 *
 * Divides accumulated densities by cell volume
 * Run after density accumulation
 */
kernel void normalize_density(
    device atomic_float* baryon_density [[buffer(0)]],
    device atomic_float* charge_density [[buffer(1)]],
    constant DensityParams& params [[buffer(2)]],
    uint cell_idx [[thread_position_in_grid]])
{
    uint n_cells = params.nx * params.ny * params.nz;
    if (cell_idx >= n_cells) return;

    float cell_volume = params.dx * params.dx * params.dx;
    float inv_volume = 1.0f / cell_volume;

    /* Load, normalize, store */
    float baryon = atomic_load_explicit(&baryon_density[cell_idx],
                                       memory_order_relaxed);
    float charge = atomic_load_explicit(&charge_density[cell_idx],
                                       memory_order_relaxed);

    baryon *= inv_volume;
    charge *= inv_volume;

    atomic_store_explicit(&baryon_density[cell_idx], baryon,
                         memory_order_relaxed);
    atomic_store_explicit(&charge_density[cell_idx], charge,
                         memory_order_relaxed);
}
