#ifndef _ASU_UTILS_H_
#define _ASU_UTILS_H_
#include <stdio.h>
#include "asu.h"
#include "molecule.h"

/*
Allocates an asymmetric unit for the given molecule types and
stoichiometry, fills the species and van der Waals radii and builds the
molecule index tables. Coordinates are left unset.

Args:
    unit: Uninitialized asymmetric unit; freed with asu_free() on success.
    mol: One molecule per type (length n_mol_types).
    stoic: Copies of each molecule type (length n_mol_types, all >= 1,
        summing to at least 2: a single molecule has no pair to score).
    n_mol_types: Number of molecule types (>= 1).

Returns:
    0 on success, -1 on invalid input (message printed to stderr).
*/
int asu_init(asu *unit, const molecule *mol, const int *stoic, int n_mol_types);

// Releases everything allocated by asu_init(). Safe to call once only.
void asu_free(asu *unit);

// Shifts the asymmetric unit so that its geometric centre is at the origin.
void asu_recenter(asu *unit);

// Prints the atoms in geometry.in style to stdout (debugging aid).
void asu_print(const asu *unit);

/*
Writes one asymmetric unit as a geometry.out block, in the same layout as
the crystal generator: "####### BEGIN STRUCTURE #######", "#key = value"
metadata lines, one "atom x y z species" line per atom, and
"#######  END  STRUCTURE #######". No lattice_vector lines are written.
Metadata: structure_number, number_of_atoms, number_of_molecules,
number_of_molecule_types, stoichiometry, number_of_atoms_in_molecule_type,
molecule_types, molecule_index (first atom of each molecule) and sr.

Args:
    unit: The asymmetric unit to serialize.
    out: Open, writable file handle (blocks are appended in order).
    structure_number: 1-based index of this structure within the file.
*/
void asu_write_block(const asu *unit, FILE *out, int structure_number);

#endif
