#ifndef _CHECK_STRUCTURE_H
#define _CHECK_STRUCTURE_H

#include "molecule.h"
#include "crystal.h"

typedef struct
{
	float (*L)[3];
	float *com1;
	float *com2;
	float *X;
	float *Y;
	float *Z;
	int index1;
	int num_atoms1;
	int index2;
	int num_atoms2;

}xtal_molecule_pair;


int check_structure(crystal random_crystal, float sr);
int check_structure_with_vdw_matrix(crystal random_crystal,
	float *vdw_matrix,
	int dim1,
	int dim2);

//void vector_cpy(float A[], float B[][3], int index);
float pdist(float T[3][3],
			float T_inv[3][3],
			float x1,
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  );

int structure_checker(crystal *xtal,
	float *vdw_cutoff,
	int total_atoms,
	int *mol_id,
	int num_mols
	);

float find_mol_len(float *X, float *Y, float *Z, int len);

void find_mol_com(float *X, float *Y, float *Z, int len, float com[3]);

float pdist_appx(float T[3][3],
			float T_inv[3][3],
			float x1,
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  );

void create_vdw_matrix_from_sr( molecule *mol,
								float *vdw_matrix,
								float sr,
								int Z);

// Looks up the van der Waals radius of a species given as two
// chars, space padded ("C ", "Cl"). Returns 0 and sets *radius, or -1 if
// the species is not in the table. Never exits; use at API boundaries.
int atom_vdw_radius(char c0, char c1, float *radius);

// Fills atom_vdw with the van der Waals radius of each atom,
// where atom holds two chars per atom. Exits on an unknown species.
void convert_atom2atom_vdw(char *atom, float *atom_vdw, int num_atoms);



#endif
