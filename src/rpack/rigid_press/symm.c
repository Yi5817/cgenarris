
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "d_algebra.h"
#include "../../crystal.h"
#include "../../spglib.h"


void symmetrize_crystal(crystal *xtal);
void get_fractional_coords(crystal *xtal, double *frac);

// basic version 1 
void symmetrize_crystal(crystal *xtal)
{
    const int spg = xtal->spg;
    
    // Get symm operations
    double translations[192][3];
    int rotations[192][3][3];
    int hall_number = hall_number_from_spg(spg);
    int Z = spg_get_symmetry_from_database(rotations,
					   translations,
					   hall_number);

    // Useful variables
    const int N = xtal->num_atoms_in_molecule;
    const int n_atoms = N * Z;
    double lat_vec_trans[3][3];
    double lat_vec[3][3] ={
	{xtal->lattice_vectors[0][0], xtal->lattice_vectors[0][1], xtal->lattice_vectors[0][2]},
	{xtal->lattice_vectors[1][0], xtal->lattice_vectors[1][1], xtal->lattice_vectors[1][2]},
	{xtal->lattice_vectors[2][0], xtal->lattice_vectors[2][1], xtal->lattice_vectors[2][2]}
    };
    mat3b3_transpose(lat_vec_trans, lat_vec);

    // First convert to fractional system
    double frac[n_atoms];
    get_fractional_coords(xtal, frac);

    // Strictly enforce symmetry by generating full crystal from asymm unit
    for(int op = 0; op < Z; op++)
    {
        double trans[3] = {translations[op][0],
			  translations[op][1],
			  translations[op][2]};
        int (*rot)[3] = rotations[op];

	// Loop over atoms in asymmetric molecule and
	// generate the atoms in the mol corresponding to op
        for(int asym = 0; asym < N; asym++)
	{
	    double *atom = frac + 3*asym;
	    
	    // Apply the symm operation to atom
            vector3_intmat3b3_multiply(rot, atom, atom);
            vector3_add(trans, atom, atom);
            // Convert to cartesian
            vector3_mat3b3_multiply(lat_vec_trans, atom, atom);
	    
	    xtal->Xcord[op*N + asym] = atom[0];
	    xtal->Ycord[op*N + asym] = atom[1];
	    xtal->Zcord[op*N + asym] = atom[2];
	}
    }
}


//version 2
/*
void symmetrize_crystal_2(crystal *xtal)
{
    const double tol = 0.001;
    const int spg = xtal->spg;
    
    // Get symm operations
    double trans[192][3];
    int rot[192][3][3];
    int hall_number = hall_number_from_spg(spg);
    int Z = spg_get_symmetry_from_database(rot,
					   trans,
					   hall_number);

    // First convert to fractional system
    const int N = xtal->num_atoms_in_molecule;
    const int n_atoms = N * Z;
    double frac[n_atoms];
    get_fractional_coords(xtal, frac);

    // Find all atoms that are symmetry-related to the atoms in the
    // asymmetric molecule.
    size_t symm_related[n_atoms];

    for(int i = 0; i < n_atoms; i++)
    {
	if (i < N)
	{
	    symm_related[i] = i;
	    break;
	}

	// find the symm equivalent atom of the asymm molecule
	double i_pos[3] = {*(frac + 3*i + 0), *(frac + 3*i + 1), *(frac + 3*i + 2)};
	for(int j = 0; j < N; j++)
	{
	    double j_pos[3] = {*(frac + 3*j + 0), *(frac + 3*j + 1), *(frac + 3*j + 2)};
	    vector3_intmat3b3_multiply(rot, atom, atom);
	    vector3_add(trans, atom, atom);

	}
    }
    
}
*/
int main()
{
    return 0;

}

void get_fractional_coords(crystal *xtal, double *frac)
{
    double lattice_vectors[3][3] ={
	{xtal->lattice_vectors[0][0], xtal->lattice_vectors[0][1], xtal->lattice_vectors[0][2]},
	{xtal->lattice_vectors[1][0], xtal->lattice_vectors[1][1], xtal->lattice_vectors[1][2]},
	{xtal->lattice_vectors[2][0], xtal->lattice_vectors[2][1], xtal->lattice_vectors[2][2]}
    };
    double inv_lattice_vectors[3][3];
    inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
    mat3b3_transpose(inv_lattice_vectors, inv_lattice_vectors);

    int total_atoms = xtal->Z *xtal->num_atoms_in_molecule;

    for(int i = 0; i < total_atoms; i++)
    {
        double temp[3] = {xtal->Xcord[i], xtal->Ycord[i], xtal->Zcord[i] };
        vector3_mat3b3_multiply(inv_lattice_vectors, temp, temp );
	vector3_frac(temp);
	*(frac + 3*i + 0) = temp[0];
	*(frac + 3*i + 1) = temp[1];
	*(frac + 3*i + 2) = temp[2];
    }
}
