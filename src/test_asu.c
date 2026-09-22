/*
 * Unit and end-to-end tests for asymmetric-unit (ASU) generation.
 *
 * Run from a directory that holds geometry_0.in and geometry_1.in
 * (JESMOK). Exits 0 when every check passes.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "asu.h"
#include "asu_generation.h"
#include "asu_utils.h"
#include "molecule.h"
#include "molecule_utils.h"
#include "randomgen.h"
#include "read_input.h"

#define N_ASU        60
#define SR_MIN       0.75
#define SR_MAX       1.30
#define MAX_ATTEMPTS 100000L
#define SEED         1234

static int n_failed = 0;

static void expect(int ok, const char *what)
{
    printf("%s: %s\n", ok ? "PASS" : "FAIL", what);
    if(!ok)
        n_failed++;
}

// Builds a molecule with num_atoms atoms of one species at the given points.
static molecule make_molecule(int num_atoms, const char *species,
                              const float xyz[][3])
{
    molecule m;
    m.num_of_atoms = num_atoms;
    m.X = (float *)malloc(num_atoms * sizeof(float));
    m.Y = (float *)malloc(num_atoms * sizeof(float));
    m.Z = (float *)malloc(num_atoms * sizeof(float));
    m.atoms = (char *)malloc(2 * num_atoms * sizeof(char));
    for(int i = 0; i < num_atoms; i++)
    {
        m.X[i] = xyz[i][0];
        m.Y[i] = xyz[i][1];
        m.Z[i] = xyz[i][2];
        m.atoms[2 * i] = species[0];
        m.atoms[2 * i + 1] = species[1];
    }
    return m;
}

static void free_molecule(molecule *m)
{
    free(m->X);
    free(m->Y);
    free(m->Z);
    free(m->atoms);
}

// asu_init: block layout for stoic = {2, 1} with 3- and 2-atom molecules.
static void test_init_indexing(void)
{
    const float a[3][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    const float b[2][3] = {{0, 0, 0}, {1.4f, 0, 0}};
    molecule mol[2] = {make_molecule(3, "C ", a), make_molecule(2, "Cl", b)};
    int stoic[2] = {2, 1};

    asu unit;
    expect(asu_init(&unit, mol, stoic, 2) == 0, "asu_init accepts stoic {2,1}");
    expect(unit.n_mols == 3 && unit.n_atoms == 8, "asu_init counts 3 molecules, 8 atoms");
    expect(unit.mol_types[0] == 0 && unit.mol_types[1] == 0 && unit.mol_types[2] == 1,
           "asu_init orders molecules by type then copy");
    expect(unit.mol_index[0] == 0 && unit.mol_index[1] == 3 && unit.mol_index[2] == 6,
           "asu_init atom offsets are 0, 3, 6");
    expect(unit.atoms[2 * 6] == 'C' && unit.atoms[2 * 6 + 1] == 'l',
           "asu_init copies species per block");
    expect(fabsf(unit.vdw_radii[0] - 1.7f) < 1e-6f && fabsf(unit.vdw_radii[7] - 1.75f) < 1e-6f,
           "asu_init fills vdW radii from the species table");
    asu_free(&unit);

    int bad_stoic[2] = {1, 0};
    expect(asu_init(&unit, mol, bad_stoic, 2) == -1, "asu_init rejects stoichiometry 0");

    molecule unknown = make_molecule(1, "Xx", a);
    expect(asu_init(&unit, &unknown, stoic, 1) == -1, "asu_init rejects an unknown species");
    free_molecule(&unknown);

    free_molecule(&mol[0]);
    free_molecule(&mol[1]);
}

static void test_min_sr(void)
{
    const float c[1][3] = {{0, 0, 0}};
    molecule mol[1] = {make_molecule(1, "C ", c)};
    int stoic[1] = {2};
    asu unit;
    asu_init(&unit, mol, stoic, 1);

    unit.Xcord[0] = 0;   unit.Ycord[0] = 0; unit.Zcord[0] = 0;
    unit.Xcord[1] = 3.4f; unit.Ycord[1] = 0; unit.Zcord[1] = 0;
    expect(fabsf(asu_min_sr(&unit) - 1.0f) < 1e-5f, "asu_min_sr: C..C at 3.4 A gives sr = 1");

    unit.Xcord[1] = 1.7f;
    expect(fabsf(asu_min_sr(&unit) - 0.5f) < 1e-5f, "asu_min_sr: C..C at 1.7 A gives sr = 0.5");

    asu_free(&unit);
    free_molecule(&mol[0]);
}

// A single-molecule "unit" has no pair to score and is rejected at init.
static void test_single_molecule_rejected(void)
{
    const float a[3][3] = {{-1, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    molecule mol[1] = {make_molecule(3, "C ", a)};
    int stoic[1] = {1};
    asu unit;
    expect(asu_init(&unit, mol, stoic, 1) == -1,
           "asu_init rejects a single molecule (sum of stoichiometry < 2)");
    free_molecule(&mol[0]);
}

// Box length: a = (3 * sum stoic_i d_i^3)^(1/3) with d_i = 2 * max |r - com|.
static void test_box_length(void)
{
    const float a[2][3] = {{-1, 0, 0}, {1, 0, 0}};      // d = 2
    const float b[2][3] = {{-2, 0, 0}, {2, 0, 0}};      // d = 4
    molecule mol[2] = {make_molecule(2, "C ", a), make_molecule(2, "C ", b)};
    int stoic[2] = {2, 1};
    float expected = cbrtf(3.0f * (2 * 8 + 1 * 64));
    expect(fabsf(asu_box_length(mol, stoic, 2) - expected) < 1e-4f,
           "asu_box_length matches (3 sum n_i d_i^3)^(1/3)");
    free_molecule(&mol[0]);
    free_molecule(&mol[1]);
}

// Every accepted unit has sr inside the window and is centred at the origin.
static void test_generate_one_window(void)
{
    molecule mol[2];
    read_geometry(&mol[0], "geometry_0.in");
    read_geometry(&mol[1], "geometry_1.in");
    recenter_molecule(&mol[0]);
    recenter_molecule(&mol[1]);
    int stoic[2] = {1, 1};
    asu unit;
    asu_init(&unit, mol, stoic, 2);
    float box_len = asu_box_length(mol, stoic, 2);
    init_genrand(SEED);

    int ok_window = 1, ok_centre = 1, ok_attempts = 1;
    for(int n = 0; n < 20; n++)
    {
        long used = asu_generate_one(&unit, mol, box_len, SR_MIN, SR_MAX, MAX_ATTEMPTS);
        if(used < 1)
            ok_attempts = 0;
        float sr = asu_min_sr(&unit);
        if(!(sr > SR_MIN && sr < SR_MAX) || fabsf(sr - unit.sr) > 1e-5f)
            ok_window = 0;
        float com[3] = {0, 0, 0};
        for(int i = 0; i < unit.n_atoms; i++)
        {
            com[0] += unit.Xcord[i];
            com[1] += unit.Ycord[i];
            com[2] += unit.Zcord[i];
        }
        if(fabsf(com[0]) + fabsf(com[1]) + fabsf(com[2]) > 1e-3f * unit.n_atoms)
            ok_centre = 0;
    }
    expect(ok_attempts, "asu_generate_one finds a unit within the attempt budget");
    expect(ok_window, "asu_generate_one: sr inside the window and stored in unit.sr");
    expect(ok_centre, "asu_generate_one recentres the unit at the origin");

    // An impossible window is reported as 0 attempts used.
    expect(asu_generate_one(&unit, mol, box_len, 5.0f, 5.1f, 50) == 0,
           "asu_generate_one returns 0 when the window is unattainable");

    asu_free(&unit);
    free_molecule(&mol[0]);
    free_molecule(&mol[1]);
}

/*
Counts geometry.out blocks in path. Every block must have n_atoms atom
lines, sr inside (sr_min, sr_max), a molecule_index line and consecutive
structure numbers starting at 1. Returns the block count, or -1.
*/
static int check_blocks(const char *path, int n_atoms, float sr_min, float sr_max)
{
    FILE *f = fopen(path, "r");
    if(!f)
        return -1;
    char line[4096];
    int n_blocks = 0, bad = 0;
    int atoms_seen = 0, number = 0, has_index = 0;
    float sr = -1;
    while(fgets(line, sizeof(line), f))
    {
        if(!strncmp(line, "####### BEGIN", 13))
        {
            atoms_seen = 0; number = 0; has_index = 0; sr = -1;
        }
        else if(!strncmp(line, "#structure_number = ", 20))
            number = atoi(line + 20);
        else if(!strncmp(line, "#molecule_index = 0 ", 20))
            has_index = 1;
        else if(!strncmp(line, "#sr = ", 6))
            sr = (float)atof(line + 6);
        else if(!strncmp(line, "atom ", 5))
            atoms_seen++;
        else if(!strncmp(line, "#######  END", 12))
        {
            n_blocks++;
            if(atoms_seen != n_atoms || number != n_blocks || !has_index ||
               !(sr > sr_min && sr < sr_max))
                bad++;
        }
    }
    fclose(f);
    return bad ? -1 : n_blocks;
}

// MPI driver: all requested units are written and pass the window check.
static void test_generate_mpi(MPI_Comm comm)
{
    molecule mol[2];
    read_geometry(&mol[0], "geometry_0.in");
    read_geometry(&mol[1], "geometry_1.in");
    int stoic[2] = {1, 1};
    int n_atoms = mol[0].num_of_atoms + mol[1].num_of_atoms;

    int n = asu_generate_mpi(mol, 2, stoic, N_ASU, SR_MIN, SR_MAX, MAX_ATTEMPTS,
                             SEED, "asu_test.out", comm);
    expect(n == N_ASU, "asu_generate_mpi generates the requested number of units");
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);
    if(my_rank == 0)
        expect(check_blocks("asu_test.out", n_atoms, SR_MIN, SR_MAX) == N_ASU,
               "asu_generate_mpi writes every unit with the right atoms, sr and numbering");

    expect(asu_generate_mpi(mol, 2, stoic, N_ASU, 1.3, 0.75, MAX_ATTEMPTS,
                            SEED, "asu_test.out", comm) == -1,
           "asu_generate_mpi rejects sr_min >= sr_max");
    expect(asu_generate_mpi(mol, 2, stoic, N_ASU, SR_MIN, SR_MAX, 0,
                            SEED, "asu_test.out", comm) == -1,
           "asu_generate_mpi rejects max_attempts < 1");
    expect(asu_generate_mpi(mol, 2, stoic, 0, SR_MIN, SR_MAX, MAX_ATTEMPTS,
                            SEED, "asu_test.out", comm) == 0,
           "asu_generate_mpi with num_structures 0 writes nothing and returns 0");

    free_molecule(&mol[0]);
    free_molecule(&mol[1]);
}

// Flat-array front end: shape and species-length errors return -1.
static void test_flat_interface(MPI_Comm comm)
{
    double positions[3 * 3] = {0, 0, 0, 1.2, 0, 0, 0, 0, 0};
    int n_atoms_per_mol[2] = {2, 1};
    int stoic[2] = {1, 1};

    expect(generate_asymmetric_units(positions, 3, 3, "C H C ", n_atoms_per_mol, 2,
                                     stoic, 2, 4, SR_MIN, SR_MAX, MAX_ATTEMPTS, SEED,
                                     "asu_flat.out", comm) == 4,
           "generate_asymmetric_units: valid flat input generates 4 units");
    expect(generate_asymmetric_units(positions, 3, 2, "C H C ", n_atoms_per_mol, 2,
                                     stoic, 2, 4, SR_MIN, SR_MAX, MAX_ATTEMPTS, SEED,
                                     "asu_flat.out", comm) == -1,
           "generate_asymmetric_units rejects positions without 3 columns");
    expect(generate_asymmetric_units(positions, 3, 3, "C H C", n_atoms_per_mol, 2,
                                     stoic, 2, 4, SR_MIN, SR_MAX, MAX_ATTEMPTS, SEED,
                                     "asu_flat.out", comm) == -1,
           "generate_asymmetric_units rejects a species string of the wrong length");
    expect(generate_asymmetric_units(positions, 3, 3, "C H C ", n_atoms_per_mol, 2,
                                     stoic, 1, 4, SR_MIN, SR_MAX, MAX_ATTEMPTS, SEED,
                                     "asu_flat.out", comm) == -1,
           "generate_asymmetric_units rejects mismatched stoichiometry length");
    int one_copy[2] = {1, 0};
    expect(generate_asymmetric_units(positions, 3, 3, "C H C ", n_atoms_per_mol, 2,
                                     one_copy, 2, 4, SR_MIN, SR_MAX, MAX_ATTEMPTS, SEED,
                                     "asu_flat.out", comm) == -1,
           "generate_asymmetric_units rejects a zero stoichiometry entry");
    int wrong_atoms[2] = {2, 2};
    expect(generate_asymmetric_units(positions, 3, 3, "C H C ", wrong_atoms, 2,
                                     stoic, 2, 4, SR_MIN, SR_MAX, MAX_ATTEMPTS, SEED,
                                     "asu_flat.out", comm) == -1,
           "generate_asymmetric_units rejects n_atoms_per_mol not summing to n_atoms");
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    if(my_rank == 0)
    {
        test_init_indexing();
        test_min_sr();
        test_single_molecule_rejected();
        test_box_length();
        test_generate_one_window();
    }
    MPI_Barrier(comm);
    test_generate_mpi(comm);
    test_flat_interface(comm);

    int total_failed = 0;
    MPI_Allreduce(&n_failed, &total_failed, 1, MPI_INT, MPI_SUM, comm);
    if(my_rank == 0)
        printf(total_failed ? "ASU tests: %d check(s) failed\n"
                            : "ASU tests: all checks passed\n", total_failed);
    MPI_Finalize();
    return total_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
