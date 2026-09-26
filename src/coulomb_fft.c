#include "grid.h"
#include <math.h>

#ifdef USE_FFTW
#include <fftw3.h>
#endif

/* FFT-based Poisson solver for Coulomb potential */
void grid_solve_coulomb_fft(Grid *grid) {
#ifdef USE_FFTW
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    float dx = grid->dx;

    /* Allocate FFT arrays */
    double *rho_in = fftw_alloc_real(nx * ny * nz);
    fftw_complex *rho_out = fftw_alloc_complex(nx * ny * (nz/2 + 1));

    /* Copy charge density to FFT input */
    for (int i = 0; i < nx * ny * nz; i++) {
        rho_in[i] = grid->charge_density[i];
    }

    /* Forward FFT */
    fftw_plan forward = fftw_plan_dft_r2c_3d(nx, ny, nz, rho_in, rho_out, FFTW_ESTIMATE);
    fftw_execute(forward);

    /* Solve in k-space: φ(k) = 4πα * ρ(k) / k² */
    float alpha = 1.44f;  /* e² in MeV·fm */
    float factor = 4.0f * PI * alpha * dx * dx * dx;

    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int iz = 0; iz < nz/2 + 1; iz++) {
                /* Wave vector */
                float kx = (ix <= nx/2) ? ix : ix - nx;
                float ky = (iy <= ny/2) ? iy : iy - ny;
                float kz = iz;

                kx *= 2.0f * PI / (nx * dx);
                ky *= 2.0f * PI / (ny * dx);
                kz *= 2.0f * PI / (nz * dx);

                float k2 = kx*kx + ky*ky + kz*kz;

                int idx = (ix * ny + iy) * (nz/2 + 1) + iz;

                if (k2 > 1e-10f) {
                    rho_out[idx][0] *= factor / k2;
                    rho_out[idx][1] *= factor / k2;
                } else {
                    /* k=0 component (constant offset) - set to zero */
                    rho_out[idx][0] = 0.0;
                    rho_out[idx][1] = 0.0;
                }
            }
        }
    }

    /* Inverse FFT */
    fftw_plan backward = fftw_plan_dft_c2r_3d(nx, ny, nz, rho_out, rho_in, FFTW_ESTIMATE);
    fftw_execute(backward);

    /* Copy result to Coulomb potential grid (normalize by N) */
    double norm = 1.0 / (nx * ny * nz);
    for (int i = 0; i < nx * ny * nz; i++) {
        grid->coulomb_potential[i] = rho_in[i] * norm;
    }

    /* Cleanup */
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_free(rho_in);
    fftw_free(rho_out);

#else
    /* Fallback: use direct summation (slow but works) */
    grid_solve_coulomb(grid);
#endif
}
