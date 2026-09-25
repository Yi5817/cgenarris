#ifndef _ASU_H_
#define _ASU_H_

/*
asu = asymmetric unit: a cluster of at least two whole, rigid molecules
with no lattice.
*/
typedef struct
{
    float *Xcord;            // n_atoms
    float *Ycord;            // n_atoms
    float *Zcord;            // n_atoms
    char  *atoms;            // 2 chars per atom, space padded: "C ", "Cl"
    float *vdw_radii;        // n_atoms, from the species table
    int   *mol_index;        // n_mols: first atom of each molecule
    int   *mol_types;        // n_mols: type of each molecule
    int   *n_atoms_in_mol;   // n_mol_types
    int   *stoic;            // n_mol_types: copies of each type
    int   n_mols;            // sum of stoic, >= 2
    int   n_mol_types;
    int   n_atoms;
    float sr;                // min d_ij / (r_i + r_j) over molecule pairs
} asu;

#endif
