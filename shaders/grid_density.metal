#include <metal_stdlib>
using namespace metal;

/* Particle structure matching CPU-side Particle */
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

/* Grid parameters passed as uniform buffer */
struct GridParams {
    int nx, ny, nz;
    float dx;
    float origin_x, origin_y, origin_z;
    int ff_range;
};

/* Form factor (Monaghan cubic spline kernel) */
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

/* Grid indexing */
static inline int grid_index(const GridParams params, int ix, int iy, int iz) {
    return ix * params.ny * params.nz + iy * params.nz + iz;
}

/* Convert world coordinates to grid cell indices */
static inline void grid_coords(const GridParams params, float x, float y, float z,
                                thread int* ix, thread int* iy, thread int* iz) {
    *ix = (int)((x - params.origin_x) / params.dx + 0.5f);
    *iy = (int)((y - params.origin_y) / params.dx + 0.5f);
    *iz = (int)((z - params.origin_z) / params.dx + 0.5f);
}

/*
 * Grid density calculation kernel
 *
 * Threading: One thread per particle
 * Each thread:
 *   1. Finds particle's grid cell
 *   2. Loops over ±ff_range neighborhood (typically ±2 cells)
 *   3. Evaluates separable form factor
 *   4. Atomically accumulates to grid density arrays
 *
 * Memory access:
 *   - Reads: particles (coalesced per threadgroup)
 *   - Writes: baryon_density, charge_density (atomic, scattered)
 */
kernel void calculate_grid_density(
    device const Particle* particles [[buffer(0)]],
    device atomic_float* baryon_density [[buffer(1)]],
    device atomic_float* charge_density [[buffer(2)]],
    constant GridParams& params [[buffer(3)]],
    uint particle_idx [[thread_position_in_grid]],
    uint num_particles [[threads_per_grid]])
{
    /* Bounds check */
    if (particle_idx >= num_particles) return;

    const device Particle& p = particles[particle_idx];

    /* Find particle's home grid cell */
    int ix, iy, iz;
    grid_coords(params, p.x, p.y, p.z, &ix, &iy, &iz);

    /* Distribute over neighboring cells within form factor support */
    for (int dx = -params.ff_range; dx <= params.ff_range; dx++) {
        float ff_x = form_factor_1d((float)dx);
        if (ff_x == 0.0f) continue;

        for (int dy = -params.ff_range; dy <= params.ff_range; dy++) {
            float ff_y = form_factor_1d((float)dy);
            if (ff_y == 0.0f) continue;

            for (int dz = -params.ff_range; dz <= params.ff_range; dz++) {
                float ff_z = form_factor_1d((float)dz);
                if (ff_z == 0.0f) continue;

                int cx = ix + dx;
                int cy = iy + dy;
                int cz = iz + dz;

                /* Boundary check (no periodic BC) */
                if (cx < 0 || cx >= params.nx) continue;
                if (cy < 0 || cy >= params.ny) continue;
                if (cz < 0 || cz >= params.nz) continue;

                /* Separable form factor */
                float ff = ff_x * ff_y * ff_z;

                int idx = grid_index(params, cx, cy, cz);

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
 * Grid normalization kernel
 *
 * Threading: One thread per grid cell
 * Divides densities by cell volume
 */
kernel void normalize_grid_density(
    device atomic_float* baryon_density [[buffer(0)]],
    device atomic_float* charge_density [[buffer(1)]],
    constant GridParams& params [[buffer(2)]],
    uint cell_idx [[thread_position_in_grid]])
{
    uint n_cells = params.nx * params.ny * params.nz;
    if (cell_idx >= n_cells) return;

    float cell_volume = params.dx * params.dx * params.dx;
    float inv_volume = 1.0f / cell_volume;

    /* Load current values atomically */
    float baryon = atomic_load_explicit(&baryon_density[cell_idx],
                                       memory_order_relaxed);
    float charge = atomic_load_explicit(&charge_density[cell_idx],
                                       memory_order_relaxed);

    /* Normalize */
    baryon *= inv_volume;
    charge *= inv_volume;

    /* Store back atomically */
    atomic_store_explicit(&baryon_density[cell_idx], baryon,
                         memory_order_relaxed);
    atomic_store_explicit(&charge_density[cell_idx], charge,
                         memory_order_relaxed);
}

/*
 * Grid zeroing kernel
 *
 * Threading: One thread per grid cell
 * Zeros density arrays before accumulation
 */
kernel void zero_grid_density(
    device atomic_float* baryon_density [[buffer(0)]],
    device atomic_float* charge_density [[buffer(1)]],
    uint cell_idx [[thread_position_in_grid]],
    uint n_cells [[threads_per_grid]])
{
    if (cell_idx >= n_cells) return;

    atomic_store_explicit(&baryon_density[cell_idx], 0.0f,
                         memory_order_relaxed);
    atomic_store_explicit(&charge_density[cell_idx], 0.0f,
                         memory_order_relaxed);
}
