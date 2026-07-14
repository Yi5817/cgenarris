#ifndef _ASU_GENERATION_H_
#define _ASU_GENERATION_H_
#include "mpi.h"
#include "molecule.h"
#include "asu.h"

// Internal placement / scoring helpers.
int asu_try_one_random_place(asu *unit, float box[3][3], molecule *mol);
int generate_one_asu(asu *unit, molecule *mol, float sr_min, float sr_max, int max_attempts);
float cal_asu_sr(asu *unit);
void estimate_box(float box[3][3], int *stoic, molecule *mol, int mol_types);

/*
Generates random molecular asymmetric units in parallel over an MPI
communicator and writes them to a single extended-XYZ file.

Args:
    positions: (n_atoms_total x 3) Cartesian coordinates, row-major.
    n_atoms_total: Number of atoms across all molecule types.
    ncols: Must be 3.
    species: Two chars per atom (element symbol, space-padded).
    n_atoms_per_mol: Atom count of each molecule type (length n_mol_types).
    n_mol_types: Number of molecule types.
    stoic: Copies of each molecule type per asymmetric unit.
    n_mol_types_b: Length of stoic; must equal n_mol_types.
    num_structures: Total number of asymmetric units to generate.
    sr_min: Lower bound on the closest interatomic distance, as a fraction
        of the sum of the two atoms' van der Waals radii.
    sr_max: Upper bound on the closest interatomic distance, as a fraction
        of the sum of the two atoms' van der Waals radii.
    max_attempts: Random-placement attempt budget per rank.
    random_seed: Base RNG seed; rank r uses random_seed + r.
    output_file: Path of the merged extended-XYZ output file.
    world_comm: MPI communicator to parallelize over.

Returns:
    Total number of asymmetric units generated across all ranks.
*/
int generate_asymmetric_units(
    double *positions,
    int n_atoms_total,
    int ncols,
    char *species,
    int *n_atoms_per_mol,
    int n_mol_types,
    int *stoic,
    int n_mol_types_b,
    int num_structures,
    double sr_min,
    double sr_max,
    long max_attempts,
    int random_seed,
    char *output_file,
    MPI_Comm world_comm);

#endif
