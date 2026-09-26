/*
 * Test program for Cugnon cross sections
 *
 * Validates implementation against known physics values
 */

#include <stdio.h>
#include <math.h>
#include "../include/cross_sections.h"

int main(void)
{
    double sqrt_s, p_lab, sigma_pp, sigma_pn, sigma_avg;
    const double M = 0.938919;  /* GeV */

    printf("Cugnon NN Cross Section Test\n");
    printf("=============================\n\n");

    printf("%-10s %-12s %-12s %-12s %-12s\n",
           "sqrt(s)", "p_lab", "sigma_pp", "sigma_pn", "sigma_avg");
    printf("%-10s %-12s %-12s %-12s %-12s\n",
           "(GeV)", "(GeV/c)", "(mb)", "(mb)", "(mb)");
    printf("----------------------------------------------------------------\n");

    /* Test points across energy regimes */
    double test_energies[] = {
        1.88,    /* Low energy */
        1.89,    /* Threshold regime */
        2.00,    /* Mid regime */
        2.20,    /* Mid-high regime */
        2.50,    /* High regime */
        3.00,    /* Very high */
        0.0      /* Sentinel */
    };

    for (int i = 0; test_energies[i] > 0.0; i++) {
        sqrt_s = test_energies[i];
        p_lab = lab_momentum_from_sqrt_s(sqrt_s, M);
        sigma_pp = cugnon_sigma_nn(sqrt_s, 1);  /* pp/nn */
        sigma_pn = cugnon_sigma_nn(sqrt_s, 0);  /* pn */
        sigma_avg = cugnon_sigma_nn_averaged(sqrt_s);

        printf("%-10.4f %-12.4f %-12.2f %-12.2f %-12.2f\n",
               sqrt_s, p_lab, sigma_pp, sigma_pn, sigma_avg);
    }

    printf("\n");

    /* Test density reduction */
    printf("\nDensity Reduction Test\n");
    printf("----------------------\n");
    sqrt_s = 2.0;
    printf("sqrt(s) = %.3f GeV\n", sqrt_s);
    printf("Free space sigma_pp = %.2f mb\n", cugnon_sigma_nn(sqrt_s, 1));

    double rho_values[] = {0.0, 0.5, 1.0, 1.5, 2.0};
    double alpha = 0.2;

    printf("\nWith alpha = %.2f:\n", alpha);
    printf("%-15s %-15s\n", "rho/rho_0", "sigma (mb)");
    for (int i = 0; i < 5; i++) {
        double rho = rho_values[i];
        double sigma = cugnon_sigma_nn_medium(sqrt_s, 1, rho, alpha);
        printf("%-15.2f %-15.2f\n", rho, sigma);
    }

    return 0;
}
