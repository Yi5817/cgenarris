#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stddef.h>
#include "mpi.h"
#include "cgenarris/read_input.h"
#include "spg_generation.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "molecule_utils.h"
#include "lattice_generator.h"
#include "randomgen.h"
#include "algebra.h"
#include "cgenarris_mpi.h"
#include "cgenarris/pygenarris_mpi.h"
#include "cgenarris/asu_generation.h"

//maximum mulipicity possible
#define ZMAX 192
#define GRAIN_SIZE 10000

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

    // Every rank parses control.in so the error is reported once per rank
    // and all ranks agree on the outcome.
    Settings set;
    if(read_settings(&set, "control.in"))
    {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    if(set.generation_type == ASU)
    {
        // Molecules come from geometry_<i>.in.
        molecule *mol = (molecule *)malloc(set.n_mol_types * sizeof(molecule));
        read_molecules(mol, set.n_mol_types);
        if(my_rank == 0)
            print_input_geometries(mol, set.n_mol_types);

        int n = asu_generate_mpi(mol, set.n_mol_types, set.stoic,
                                 set.num_structures, set.asu_sr_min,
                                 set.asu_sr_max, set.max_attempts,
                                 set.random_seed, set.asu_output_file,
                                 world_comm);

        for(int i = 0; i < set.n_mol_types; i++)
        {
            free(mol[i].atoms);
            free(mol[i].X);
            free(mol[i].Y);
            free(mol[i].Z);
        }
        free(mol);
        MPI_Finalize();
        return n < 0 ? EXIT_FAILURE : EXIT_SUCCESS;
    }

    molecule *mol = (molecule*)malloc(sizeof(molecule));//store molecule
    read_geometry(mol, "geometry.in");				//read molecule from geometry.in

    int num_atoms_in_molecule = mol->num_of_atoms;
    int dim_vdw_matrix = num_atoms_in_molecule * set.Z ;
    float *vdw_cutoff_matrix = (float *) malloc( dim_vdw_matrix *
								dim_vdw_matrix *
								sizeof(float) ); //square matrix

    create_vdw_matrix_from_sr(mol, vdw_cutoff_matrix, set.sr, set.Z);


    int status = 0;
    if(set.generation_type == CRYSTAL)
    {
        status = mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
        vdw_cutoff_matrix,
        dim_vdw_matrix,
        dim_vdw_matrix,
        set.num_structures,
        set.Z,
        set.vol_mean,
        set.vol_std,
        set.tol,
        set.max_attempts,
        set.spg_dist_type,
        set.vol_attempts,
        set.random_seed,
        set.norm_dev,
        set.angle_std,
        world_comm);

    }

	else if(set.generation_type == LAYER)			// for layer generation
	{
	    mpi_generate_layer_with_vdw_cutoff_matrix(
		vdw_cutoff_matrix,
		dim_vdw_matrix,
		dim_vdw_matrix,
		set.num_structures,
		set.Z,
		set.vol_mean,
		set.vol_std,
		set.interface_area_mean,
		set.interface_area_std,
		set.volume_multiplier,
		set.tol,
		set.max_attempts,
		set.spg_dist_type,
		set.lattice_vector_2d,
		set.vol_attempts,
		set.random_seed,
		world_comm);

	}

    MPI_Finalize();
    return status < 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
