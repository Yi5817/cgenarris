#ifndef _ASU_GENERATION_H_
#define _ASU_GENERATION_H_
#include "mpi.h"
#include "molecule.h"
#include "asu.h"

/*
Random asymmetric-unit (ASU) generation.

Each molecule copy receives a random rotation about its centre and a random
translation inside a cubic box of edge asu_box_length(). A placement is
kept when the smallest intermolecular d_ij / (r_i^vdW + r_j^vdW) over all
molecule pairs, the "specific radius" sr, lies in the open window
(sr_min, sr_max). Molecules passed to these functions must already be
centred at the origin (recenter_molecule()).
*/

// Edge length a = (3 * sum_i stoic[i] * d_i^3)^(1/3) of the cubic
// placement box, where d_i is twice the largest atom-to-centre distance
// of molecule type i.
float asu_box_length(const molecule *mol, const int *stoic, int n_mol_types);

// Draws a random rotation and translation for every molecule copy and
// writes the resulting Cartesian coordinates into unit. Consumes six
// uniform random numbers per molecule copy, in molecule order.
void asu_place_random(asu *unit, const molecule *mol, float box_len);

// Smallest d_ij / (r_i + r_j) over all pairs of atoms in different
// molecules.
float asu_min_sr(const asu *unit);

/*
Repeats asu_place_random() until the sr window is met, then recentres the
unit at the origin and stores sr in unit->sr.

Returns:
    The number of placements tried (>= 1) on success, 0 if max_attempts
    placements were tried without meeting the window.
*/
long asu_generate_one(asu *unit, const molecule *mol, float box_len,
                      float sr_min, float sr_max, long max_attempts);

/*
Generates asymmetric units in parallel over an MPI communicator and writes
them to a single file of geometry.out blocks (see asu_write_block()).
Collective: every rank must call it with identical arguments.

Args:
    mol: One molecule per type (length n_mol_types). Recentred in place.
    n_mol_types: Number of molecule types (>= 1).
    stoichiometry: Copies of each molecule type per unit (all >= 1, sum >= 2).
    num_structures: Total number of asymmetric units to generate (>= 0),
        split evenly across ranks.
    sr_min: Lower bound of the sr window; 0 <= sr_min < sr_max.
    sr_max: Upper bound of the sr window.
    max_attempts: Placements tried per asymmetric unit before a rank gives
        up (>= 1). A rank that gives up stops early with fewer units.
    random_seed: Base seed; rank r seeds its generator with seed + r.
        0 draws a time-based seed on rank 0 and shares it with all ranks.
    output_file: Path of the merged output file; overwritten.
    comm: MPI communicator to parallelize over.

Returns:
    Total number of asymmetric units written across all ranks, or -1 on
    invalid input (message printed to stderr).
*/
int asu_generate_mpi(molecule *mol, int n_mol_types, const int *stoichiometry,
                     int num_structures, double sr_min, double sr_max,
                     long max_attempts, int random_seed,
                     const char *output_file, MPI_Comm comm);

/*
Flat-array front end of asu_generate_mpi() for the Python bindings.

Args:
    positions: (n_atoms_total x 3) Cartesian coordinates in Angstrom of all
        molecule types concatenated, row-major.
    n_atoms_total: Number of rows in positions; must equal the sum of
        n_atoms_per_mol.
    ncols: Number of columns in positions; must be 3.
    species: Element symbol of every atom as two chars, space padded, in
        the same order as positions ("C H Cl" -> "C H Cl"). Its length
        must be 2 * n_atoms_total.
    n_atoms_per_mol: Atom count of each molecule type (length n_mol_types).
    n_mol_types: Number of molecule types.
    stoichiometry: Copies of each molecule type per unit (all >= 1, sum >= 2).
    n_stoichiometry: Length of stoichiometry; must equal n_mol_types.
    num_structures, sr_min, sr_max, max_attempts, random_seed,
    output_file, world_comm: as in asu_generate_mpi().

Returns:
    Total number of asymmetric units written, or -1 on invalid input.
*/
int generate_asymmetric_units(
    double *positions,
    int n_atoms_total,
    int ncols,
    char *species,
    int *n_atoms_per_mol,
    int n_mol_types,
    int *stoichiometry,
    int n_stoichiometry,
    int num_structures,
    double sr_min,
    double sr_max,
    long max_attempts,
    int random_seed,
    char *output_file,
    MPI_Comm world_comm);

#endif
