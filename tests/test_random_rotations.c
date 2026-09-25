#include <math.h>
#include <stdio.h>
#include <string.h>

#include "algebra.h"
#include "randomgen.h"

#define SAMPLES 100000
#define BINS 10

static int failures;

static void expect(int ok, const char *message)
{
    if(!ok)
    {
        fprintf(stderr, "FAIL: %s\n", message);
        failures++;
    }
}

static void check_distribution(unsigned int seed_value)
{
    double mean[3][3] = {{0}}, second[3][3] = {{0}};
    int angles[BINS] = {0};
    int orientations[8*8*8] = {0};
    double max_orthogonal_error = 0, max_det_error = 0;
    const double pi = acos(-1.0);
    init_genrand(seed_value);
    for(int sample = 0; sample < SAMPLES; sample++)
    {
        float r[3][3];
        generate_random_rotation_matrix(r);
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            expect(isfinite(r[i][j]), "finite rotation entries");
            mean[i][j] += r[i][j];
            second[i][j] += (double)r[i][j]*r[i][j];
            double dot = 0;
            for(int k = 0; k < 3; k++) dot += (double)r[k][i]*r[k][j];
            max_orthogonal_error = fmax(max_orthogonal_error,
                                       fabs(dot - (i == j)));
        }
        max_det_error = fmax(max_det_error, fabs(det_mat3b3(r) - 1));
        double cosine = ((double)r[0][0] + r[1][1] + r[2][2] - 1)/2;
        double angle = acos(fmax(-1, fmin(1, cosine)));
        int bin = (int)(BINS*angle/pi);
        angles[bin < BINS ? bin : BINS-1]++;
        // In ZYZ coordinates, Haar measure is uniform in alpha, cos(beta),
        // gamma jointly. This also checks twist, unlike a direction-only test.
        double coordinates[3] = {
            (atan2(r[1][2], r[0][2]) + pi)/(2*pi),
            ((double)r[2][2] + 1)/2,
            (atan2(r[2][1], -r[2][0]) + pi)/(2*pi)
        };
        int cell = 0;
        for(int j = 0; j < 3; j++)
        {
            int index = (int)(8*fmax(0, fmin(1, coordinates[j])));
            cell = 8*cell + (index < 8 ? index : 7);
        }
        orientations[cell]++;
    }
    expect(max_orthogonal_error < 2e-6, "R transpose R equals identity");
    expect(max_det_error < 2e-6, "determinant equals +1");
    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
        // Every rotated coordinate axis is uniform on the unit sphere.
        expect(fabs(mean[i][j]/SAMPLES) < .006, "zero orientation mean");
        expect(fabs(second[i][j]/SAMPLES - 1.0/3) < .006,
               "isotropic second moment");
    }
    int cumulative = 0;
    for(int bin = 0; bin < BINS-1; bin++)
    {
        cumulative += angles[bin];
        double theta = (bin+1)*pi/BINS;
        // Haar rotation-angle CDF on [0, pi]: (theta - sin(theta))/pi.
        expect(fabs((double)cumulative/SAMPLES - (theta-sin(theta))/pi) < .006,
               "Haar rotation-angle distribution");
    }
    double expected = (double)SAMPLES/512, chi_square = 0;
    for(int i = 0; i < 512; i++)
    {
        double delta = orientations[i] - expected;
        chi_square += delta*delta/expected;
    }
    // chi-square(511) upper-tail critical value at 1e-6 is 677.59972.
    expect(chi_square < 677.6, "joint uniformity of complete orientations");
    printf("seed %u: joint-orientation chi-square = %.3f (511 df)\n",
           seed_value, chi_square);
    printf("seed %u: mean R00 = %.6f, max orthogonality error = %.3g\n",
           seed_value, mean[0][0]/SAMPLES, max_orthogonal_error);
}

int main(void)
{
    check_distribution(1234);
    check_distribution(5678);
    float first[3][3], repeat[3][3];
    init_genrand(42);
    generate_random_rotation_matrix(first);
    float next = uniform_dist_01();
    init_genrand(42);
    generate_random_rotation_matrix(repeat);
    expect(memcmp(first, repeat, sizeof(first)) == 0, "seed reproducibility");
    init_genrand(42);
    for(int i = 0; i < 3; i++) (void)uniform_dist_01();
    expect(next == uniform_dist_01(), "exactly three RNG draws per rotation");
    return failures ? 1 : 0;
}
