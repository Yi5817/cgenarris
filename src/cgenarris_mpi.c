#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stddef.h>
#include "mpi.h"
#include "read_input.h"
#include "spg_generation.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "molecule_utils.h"
#include "lattice_generator.h"
#include "randomgen.h"
#include "algebra.h"
#include "cgenarris_mpi.h"
#include "pygenarris_mpi.h"

//maximum mulipicity possible
#define ZMAX 192
#define GRAIN_SIZE 10000

extern int *seed;
extern unsigned int *seed2;
extern float TOL;

void create_vdw_matrix_from_sr(molecule *mol,
			       float *vdw_matrix,
			       float sr,
			       int Z);

int main(int argc, char **argv)
{
	//Initialise MPI
    MPI_Init(&argc, &argv);
    int total_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &total_ranks);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm world_comm = MPI_COMM_WORLD;

    //variable declarartion
    molecule *mol = (molecule*)malloc(sizeof(molecule));//store molecule
    //////added here
    char generation_type[25]; //indicate whether its molecular crystal generation or layer generation
    float interface_area_mean; //for layer
    float interface_area_std; //for layer
    int volume_multiplier; // for layer
    float  lattice_vector_2d[2][3]; //for layer
    float volume_std;	//standard dev for volumes
    float volume_mean;	//mean volume
    float sr;			//specific radius proportion for structure checks
    float Zp_max;		//Z'' . not implemented
    int Z;				//multiplicity of general position
    int num_structures;	//num of structures per spg
    int rigid_press; // Whether rigid press optimization should be performed
    long max_attempts;	//max attempts per space group
    float tol;
    char spg_dist_type[10];  //spg distribution type
    int vol_attempt;    // no of attempts after which volume is resampled
    int random_seed;    //seed for random gen
    float norm_dev;
    float angle_std;
    int mol_types;  // No. of distinct molecule types
    int *stoic;  // kept for compatibility with read_control()

    read_geometry(mol, "geometry.in");				//read molecule from geometry.in
    read_control(&num_structures,
       		 &Z,
       		 &Zp_max,
       		 &volume_mean,
       		 &volume_std,
       		 &sr,
       		 &max_attempts,
                 spg_dist_type,
                 &vol_attempt,
                 &random_seed,
       		 generation_type,
       		 &interface_area_mean,
       		 &interface_area_std,
       		 &volume_multiplier,
       		 lattice_vector_2d,
                 &norm_dev,
                 &angle_std,
                 &stoic,
		 &mol_types,
		 &rigid_press);	//get settings

    tol = TOL;
    int num_atoms_in_molecule = mol->num_of_atoms;
    int dim_vdw_matrix = num_atoms_in_molecule * Z ;
    float *vdw_cutoff_matrix = (float *) malloc( dim_vdw_matrix *
								dim_vdw_matrix *
								sizeof(float) ); //square matrix

    create_vdw_matrix_from_sr(mol, vdw_cutoff_matrix, sr, Z);


    if(!strcmp(generation_type, "crystal"))
    {
        mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
        vdw_cutoff_matrix,
        dim_vdw_matrix,
        dim_vdw_matrix,
        num_structures,
        Z,
        volume_mean,
        volume_std,
        tol,
        max_attempts,
        spg_dist_type,
        vol_attempt,
        random_seed,
        norm_dev,
        angle_std,
        world_comm);

    }

	else if(!strcmp(generation_type, "layer"))			// for layer generation
	{
	    mpi_generate_layer_with_vdw_cutoff_matrix(
		vdw_cutoff_matrix,
		dim_vdw_matrix,
		dim_vdw_matrix,
		num_structures,
		Z,
		volume_mean,
		volume_std,
		interface_area_mean,
		interface_area_std,
		volume_multiplier,
		tol,
		max_attempts,
		spg_dist_type,
		lattice_vector_2d,
		vol_attempt,
		random_seed,
		world_comm);

	}

    else
    {
        printf("***ERROR: generation type unsupported : %s\n", generation_type);
    }

    MPI_Finalize();
    return 0;
}
