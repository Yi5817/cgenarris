/*
 * Standalone test for lattice_generator.
 *
 * For each crystal system, generates many random lattices and checks that
 * the cell volume matches the target and that the lattice parameters obey
 * the constraints of that system. Prints PASS / FAIL; exit code 0 on PASS.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lattice_generator.h"
#include "randomgen.h"

#define N_TRIALS      2000
#define TARGET_VOLUME 1234.5f
#define NORM_STD      0.3f
#define ANGLE_STD     10.0f
#define REL_TOL       1e-3
#define ANGLE_TOL     1e-2   /* degrees */

/* defined in pygenarris*.c in the main programs; provide it here */
unsigned int *seed2;

static int n_fail = 0;

static void cell_params(float lat[3][3], double *a, double *b, double *c,
                        double *alpha, double *beta, double *gamma,
                        double *volume)
{
    double av[3] = {lat[0][0], lat[0][1], lat[0][2]};
    double bv[3] = {lat[1][0], lat[1][1], lat[1][2]};
    double cv[3] = {lat[2][0], lat[2][1], lat[2][2]};
#define DOT(u, v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])
    *a = sqrt(DOT(av, av));
    *b = sqrt(DOT(bv, bv));
    *c = sqrt(DOT(cv, cv));
    *alpha = acos(DOT(bv, cv) / (*b * *c)) * 180.0 / M_PI;
    *beta  = acos(DOT(av, cv) / (*a * *c)) * 180.0 / M_PI;
    *gamma = acos(DOT(av, bv) / (*a * *b)) * 180.0 / M_PI;
#undef DOT
    double det = av[0]*(bv[1]*cv[2] - bv[2]*cv[1])
               - av[1]*(bv[0]*cv[2] - bv[2]*cv[0])
               + av[2]*(bv[0]*cv[1] - bv[1]*cv[0]);
    *volume = fabs(det);
}

static void expect(int ok, const char *system, const char *what,
                   double got, double want)
{
    if (!ok)
    {
        n_fail++;
        if (n_fail <= 20)
            printf("FAIL %-12s %-8s got %.6f want %.6f\n",
                   system, what, got, want);
    }
}

static void check(const char *system, float lat[3][3],
                  int eq_ab, int eq_abc, int a90, int b90, int g90,
                  int g120)
{
    double a, b, c, alpha, beta, gamma, vol;
    cell_params(lat, &a, &b, &c, &alpha, &beta, &gamma, &vol);

    expect(fabs(vol - TARGET_VOLUME) / TARGET_VOLUME < REL_TOL,
           system, "volume", vol, TARGET_VOLUME);
    expect(isfinite(vol) && vol > 0, system, "finite", vol, TARGET_VOLUME);
    if (eq_ab)  expect(fabs(a - b) / a < REL_TOL, system, "a==b", a, b);
    if (eq_abc) expect(fabs(a - c) / a < REL_TOL, system, "a==c", a, c);
    if (a90) expect(fabs(alpha - 90) < ANGLE_TOL, system, "alpha", alpha, 90);
    if (b90) expect(fabs(beta - 90) < ANGLE_TOL, system, "beta", beta, 90);
    if (g90) expect(fabs(gamma - 90) < ANGLE_TOL, system, "gamma", gamma, 90);
    if (g120) expect(fabs(gamma - 120) < ANGLE_TOL, system, "gamma", gamma, 120);
    /* free angles must stay inside the generator's [30, 150] window */
    expect(alpha > 29.9 && alpha < 150.1, system, "alpha_rng", alpha, 90);
    expect(beta > 29.9 && beta < 150.1, system, "beta_rng", beta, 90);
    expect(gamma > 29.9 && gamma < 150.1, system, "gamma_rng", gamma, 90);
}

int main(void)
{
    unsigned int s = 7;
    seed2 = &s;
    init_genrand(s);

    float lat[3][3];
    for (int i = 0; i < N_TRIALS; i++)
    {
        gen_triclinic_lattice(lat, TARGET_VOLUME, NORM_STD, ANGLE_STD);
        check("triclinic", lat, 0, 0, 0, 0, 0, 0);

        gen_monoclinic_lattice(lat, TARGET_VOLUME, NORM_STD, ANGLE_STD);
        check("monoclinic", lat, 0, 0, 1, 0, 1, 0);

        gen_orthorhombic_lattice(lat, TARGET_VOLUME, NORM_STD);
        check("orthorhombic", lat, 0, 0, 1, 1, 1, 0);

        gen_tetragonal_lattice(lat, TARGET_VOLUME, NORM_STD);
        check("tetragonal", lat, 1, 0, 1, 1, 1, 0);

        gen_hexagonal_lattice(lat, TARGET_VOLUME, NORM_STD);
        check("hexagonal", lat, 1, 0, 1, 1, 0, 1);

        gen_cubic_lattice(lat, TARGET_VOLUME);
        check("cubic", lat, 1, 1, 1, 1, 1, 0);
    }

    if (n_fail)
    {
        printf("\nFAIL: %d checks failed\n", n_fail);
        return 1;
    }
    printf("PASS: all lattice checks passed (%d trials per system)\n",
           N_TRIALS);
    return 0;
}
