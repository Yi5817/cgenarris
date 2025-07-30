#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include "read_input.h"
#include "spg_generation.h"
#include "pygenarris.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "molecule_utils.h"
#include "lattice_generator.h"
#include "randomgen.h"
#include "algebra.h"
#include "spglib.h"


extern float TOL;

unsigned int *seed;
unsigned int *seed2;


/* atoms supported
 */
void create_crystal_from_array(crystal *xtal, double lattice_vector[3][3], double *Xc,int total_atoms1,
        double *Yc,int total_atoms2, double *Zc, int total_atoms3,char* atoms, int total_atoms, int Z, int spg)
{

    xtal->spg = spg;
    xtal->Z = Z;

    int num_atoms_in_molecule = total_atoms/Z;
    allocate_xtal(xtal, Z, num_atoms_in_molecule);
    xtal->num_atoms_in_molecule = num_atoms_in_molecule;

    xtal->lattice_vectors[0][0] = lattice_vector[0][0];
    xtal->lattice_vectors[0][1] = lattice_vector[0][1];
    xtal->lattice_vectors[0][2] = lattice_vector[0][2];

    xtal->lattice_vectors[1][0] = lattice_vector[1][0];
    xtal->lattice_vectors[1][1] = lattice_vector[1][1];
    xtal->lattice_vectors[1][2] = lattice_vector[1][2];

    xtal->lattice_vectors[2][0] = lattice_vector[2][0];
    xtal->lattice_vectors[2][1] = lattice_vector[2][1];
    xtal->lattice_vectors[2][2] = lattice_vector[2][2];

    for(int i = 0; i < total_atoms; i++)
    {
        xtal->Xcord[i] = Xc[i];
        xtal->Ycord[i] = Yc[i];
        xtal->Zcord[i] = Zc[i];
        xtal->atoms[2*i] = atoms[2*i];
        xtal->atoms[2*i+1] = atoms[2*i+1];

    }

    //print_crystal(xtal);
}

/*
int c_check_structure(crystal xtal, double sr)
{
    float f_sr = sr;
    return check_structure(xtal, f_sr);
}
*/


int num_compatible_spacegroups(int Z, double tolerance)
{
    //set global variable tolerance
    TOL = tolerance;

    COMPATIBLE_SPG compatible_spg[230];
    int num_compatible_spg = 0;
    int thread_num = 1;
    molecule *mol = (molecule*)malloc(sizeof(molecule));

    //read geometry from geometry.in
    read_geometry(mol);

    find_compatible_spg_positions(mol,
                                  Z,
                                  compatible_spg,
                                  &num_compatible_spg,
                                  thread_num);

    free(mol);

    return num_compatible_spg;

}
