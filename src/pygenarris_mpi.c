#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stddef.h>
#include "mpi.h"
#include "read_input.h"
#include "spg_generation.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "cocrystal_utils.h"
#include "input_settings.h"
#include "molecule_utils.h"
#include "lattice_generator.h"
#include "lattice_generator_layer.h"
#include "randomgen.h"
#include "algebra.h"
#include "pygenarris_mpi.h"
#include "pygenarris_mpi_utils.h"

#define ZMAX 192
#define BATCH_SIZE 10000
#define PI 3.141592653
#define MAX_ANGLE 150 * PI/180
#define MIN_ANGLE 30 * PI/180



unsigned int *seed;
unsigned int *seed2;
//float global_lattice_vector_2d[2][3];
//int SET_INTERFACE_AREA = 0;

extern float TOL;

void mpi_generate_cocrystals_with_vdw_matrix(
    float *vdw_matrix,
    int dim1,
    int dim2,
    int num_structures,
    int Z,
    double volume_mean1,
    double volume_std1,
    double tol1,
    long max_attempts,
    char *spg_dist_type,
    int vol_attempt,
    int random_seed,
    float norm_dev,
    float angle_std,
    int *stoic,
    int n_mol_types,
    MPI_Comm world_comm)
{

    TOL = tol1;
    float sr = -1;

    // Create settings structure
    Settings set;
    set.num_structures = num_structures;
    set.Z = Z;
    set.vol_attempts = vol_attempt;
    set.random_seed = random_seed;
    set.stoic = stoic;
    set.n_mol_types = n_mol_types;
    set.max_attempts = max_attempts;
    set.vol_mean = volume_mean1;
    set.vol_std = volume_std1;
    set.norm_dev = norm_dev;
    set.angle_std = angle_std;
    set.sr = sr;
    set.spg_dist_type = spg_dist_type;
    set.vdw_matrix = vdw_matrix;
    set.generation_type = COCRYSTAL;

    if (dim1 != dim2)
    {
        printf("***ERROR:vdw cutoff matrix is not square***\n");
        exit(EXIT_FAILURE);
    }

    //Initialise MPI
    int total_ranks;
    MPI_Comm_size(world_comm, &total_ranks);
    int my_rank;
    MPI_Comm_rank(world_comm, &my_rank);

    if(my_rank == 0)
    {
        print_input_settings(set);
        printf("Generation type :cocrystal\n");
    }

    //variable declarations
    int stop_flag = 0;  // to stop threads if limit is reached
    int success_flag = 0;  // Check if any rank generated successfully
    int struct_counter = 0;    //counts number of structures
    long attempt = 0;
    int spg;
    float volume;
    int found_poll[total_ranks];

    cocrystal *cxtal = malloc(sizeof(cocrystal));     //dummy cocrystal
    molecule *mol = malloc(n_mol_types*sizeof(molecule));//store molecules

    //file to output geometries
    FILE *out_file = open_output_file(my_rank);

    // File to output generation.log
    FILE *log_file = NULL;
    if(my_rank == 0)
        log_file = fopen("generation.log", "w");

    // Init Random seeding
    seed = malloc(sizeof(unsigned int)); //seed for uniform gen
    seed2 = malloc(sizeof(unsigned int)); //seed for random
    init_random_seed(seed, seed2, random_seed, my_rank);

    // Read molecules and recenter mol to orgin
    read_molecules(mol, n_mol_types);
    recenter_molecules(mol, n_mol_types);

    // Print molecule geometries
    if(my_rank == 0)
        print_input_geometries(mol, n_mol_types);
    MPI_Barrier(world_comm);


    // Initialzations
    int n_atoms_in_mol[n_mol_types];
    volume = draw_volume(set.vol_mean, set.vol_std);
    get_n_atoms_in_mol(n_atoms_in_mol, mol, n_mol_types);
    cxtal_init(cxtal, stoic, n_atoms_in_mol, n_mol_types, Z);

    // Find compatible space groups
    int allowed_spg[230];
    int num_allowed_spg = 0;
    find_allowed_spg(allowed_spg, &num_allowed_spg, Z);
    if(my_rank == 0)
        print_allowed_spg(allowed_spg, num_allowed_spg);

    for(int spg_index = 0; spg_index < num_allowed_spg; spg_index++)
    {
        MPI_Barrier(world_comm);

        // Get ready for generation
        spg = allowed_spg[spg_index];
        cxtal->spg = spg;
        time_t start_time = time(NULL);
        int spg_num_structures = find_num_structure_for_spg(num_structures,\
            spg_dist_type, spg, Z);  // Number of structures for spg

        if(my_rank == 0)
	    printf("Generating %d structures from spg %d...\n",
		   spg_num_structures, spg);
        //print_start_spg(spg, spg_num_structures);
        if(!spg_num_structures)
            goto end_spg_loop;

        // Loop over attempts
        // Each rank performs `BATCH_SIZE` attempts before talking to master rank
        attempt = 0;
        long max_attempt_per_rank = max_attempts/total_ranks;
        for(; 1; attempt += BATCH_SIZE)
        {
            int result = try_crystal_generation(cxtal, set,
                                                mol, &volume, attempt, BATCH_SIZE);

            // Poll to see which ranks succeded
            MPI_Gather(&result, 1, MPI_INT, found_poll, 1, MPI_INT, 0, world_comm);

            // Write structures send by slave ranks to output file
            if(my_rank == 0)
                write_structures(cxtal, found_poll, &struct_counter,
                                 spg_num_structures, out_file, total_ranks,
                                 world_comm);
            // Sends structures to master
            else
                send_structures(cxtal, result, world_comm);

            // Time to stop? -  enough structures or ran out of attempts
            stop_flag = check_stop_condition(struct_counter, spg_num_structures,
                                             attempt, max_attempt_per_rank);
            //printf("stop = %d\n", stop_flag);
            MPI_Bcast(&stop_flag, 1, MPI_INT, 0, world_comm);
            if(stop_flag)
                break;

            // Let other ranks know if a structure was generated.
            MPI_Bcast(&success_flag, 1, MPI_INT, 0, world_comm);
            if(success_flag) // Resets if a structure was generated
            {
                attempt = 0;
                volume = draw_volume(set.vol_mean, set.vol_std);  // Reset volume
            }
        }

        if(stop_flag == ATTEMPTS_STOP)
        {
            if(my_rank == 0)
	      printf("Max number of attempts %ld reached.\n", max_attempts);
            MPI_Barrier(world_comm);
        }
        else if(stop_flag == ENOUGH_STOP)
        {
            if(my_rank == 0)
                printf("Generated required number of structures..\n");
        }

        end_spg_loop:
        if(my_rank == 0)
        {
	        time_t end_time = time(NULL);
                double elapsed = difftime(end_time, start_time);
		print_spg_end(elapsed, struct_counter, spg);
        }
        struct_counter = 0;
        //break;

    }// spg generation loop

    if(my_rank == 0)
    {
        print_exit();
	fclose(out_file);
	fclose(log_file);
    }
      
}
    
void mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
    float *vdw_matrix,
    int dim1,
    int dim2,
    int num_structures,
    int Z,
    double volume_mean1,
    double volume_std1,
    double tol1,
    long max_attempts,
    char *spg_dist_type,
    int vol_attempt,
    int random_seed,
    float norm_dev,
    float angle_std,
    MPI_Comm world_comm)
{
    float Zp_max = 1;
    float volume_mean = volume_mean1;
    float volume_std = volume_std1;
    TOL = tol1;

    if (dim1 != dim2)
    {
      printf("***ERROR:vdw cutoff matrix is not square***\n");
      exit(EXIT_FAILURE);
    }

    Settings set;
    set.num_structures = num_structures;
    set.Z = Z;
    set.vol_attempts = vol_attempt;
    set.random_seed = random_seed;
    set.stoic = NULL;
    set.n_mol_types = 0;
    set.max_attempts = max_attempts;
    set.vol_mean = volume_mean1;
    set.vol_std = volume_std1;
    set.norm_dev = norm_dev;
    set.angle_std = angle_std;
    set.sr = -1;
    set.spg_dist_type = spg_dist_type;
    set.vdw_matrix = vdw_matrix;
    set.generation_type = CRYSTAL;

    //Initialise MPI
    int total_ranks;
    MPI_Comm_size(world_comm, &total_ranks);
    int my_rank;
    MPI_Comm_rank(world_comm, &my_rank);

    //variable declarations
    //int fail_count;
    int stop_flag = 0;  // to stop threads if limit is reached
    int success_flag = 0;
    int counter = 0;    //counts number of structures
    int spg_index = 0;  //space group to be generated
    FILE *out_file;     //file to output geometries
    if (my_rank == 0)
    {
        out_file = fopen("geometry.out","w");
        if(!out_file)       //check permissions
        {
            printf("***ERROR: cannot create geometry.out \n");
            exit(EXIT_FAILURE);
        }
        //fprintf(out_file, "my_rank=%d\n", my_rank);
    }
    else
    {
        out_file = NULL;
    }

    //random number seeding, different seeds for different threads
    if (random_seed == 0)
    {
        srand((unsigned int)time(NULL));
        random_seed = rand();
    }
    else
    {
        srand((unsigned int) 19023411);
    }
    seed = (unsigned int*)malloc(sizeof(unsigned int)); //seed for uniform gen
    seed2 = (unsigned int*)malloc(sizeof(unsigned int)); //seed for random
    *seed = (unsigned int)abs(my_rank*7 + random_seed);  //some random seed private for each threads
    *seed2 = (unsigned int)abs(my_rank*17 + random_seed);
    init_genrand(*seed);

    //storing information for compatible space groups
    COMPATIBLE_SPG compatible_spg[230];
    int num_compatible_spg = 0;

    //variable declarartion
    molecule *mol = (molecule*)malloc(sizeof(molecule));//store molecule
    crystal *random_crystal = (crystal*)malloc(sizeof(crystal));//dummy crystal
    //float volume_std; //standard dev for volumes
    //float volume_mean;    //mean volume
    //float sr = -1;          //specific radius proportion for structure checks
                        //see paper for definition
    //float Zp_max;     //Z'' . not implemented
    float volume;       //random volume used of generation
    //int Z;                //multiplicity of general position
    //int num_structures;   //num of structures
    int spg;            //space group attempted
    //long max_attempts;    //max attempts per space group

    //read input from file, read molecular structure from geometry,Z, Zp
    if(my_rank == 0)
    {
        int len = 100;
        char name[len];
        gethostname(name, len);
        printf("PARALLELIZATION INFO:\n");
        printf("---------------------------\n");
        printf("Using MPI for parallelization.\n");
        printf("cgenarris is running on %d processes. \n", total_ranks);
        printf("Master Rank host name: %s \n", name);
        print_time();
    }

    if(my_rank == 0)
    {
        printf("---------------------------\n");
        printf("\n");
    }
    MPI_Barrier(world_comm);

    read_geometry(mol, "geometry.in");             //read molecule from geometry.in

    //recenter molecule to origin
    recenter_molecule(mol);

    if (my_rank == 0)
    {
        print_input_geometry(mol);
	print_input_settings(set);

    }

    MPI_Barrier(world_comm);

    //inititalise volume
    do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
    int N = mol->num_of_atoms;
    allocate_xtal(random_crystal, ZMAX, N); //allcate memory
    random_crystal->num_atoms_in_molecule = mol->num_of_atoms;
    int total_atoms = mol->num_of_atoms * Z;


    if(my_rank == 0)
    {
        printf("COMPATIBLE SPACE GROUP INFO:\n");
        printf("-----------------------------------------------------\n");
        printf("Detecting compatible space groups"
               " from molecule symmetries...\n");
    }

    //compatible space groups are stored in compatible_spg. see the
    //object definition in spg_generation.h
    //every thread has its own copy
    find_compatible_spg_positions(mol,
                                  Z,
                                  compatible_spg,
                                  &num_compatible_spg,
                                  my_rank+1);

    if(my_rank == 0)
    {
        printf("\n");
        printf("Total compatible space groups = %d\n", num_compatible_spg);
        //printf("Number of molecular axes = %d \n", num_axes);
        printf("-----------------------------------------------------\n\n");
        printf("Starting generation...\n");
        print_time();
        printf("\n");
        sleep(1);
    }

    MPI_Barrier(world_comm); // wait for other friends to join

    /*deprecated
    //find allowed space groups for general position (deprecated)
    //int allowed_spg[230];
    //int num_allowed_spg = 0;
    //find_allowed_spg(allowed_spg, &num_allowed_spg, Z);
    //printf("allowed %d \n", num_compatible_spg);
    */
    while( spg_index < num_compatible_spg )
    {
        MPI_Barrier(world_comm);
        //counter counts the number of structures generated

        spg = compatible_spg[spg_index].spg; //pick a spg

        //time information
        time_t start_time = time(NULL);
        //fail_count = 0;

        // Get number of structures to be generated in this spg
        // skip to next spg if this is zero
        int spg_num_structures = find_num_structure_for_spg(num_structures, spg_dist_type, spg, Z);
        if( !spg_num_structures )
        {
            spg_index++;
            continue;
        }

        //print attempted space group
        if (my_rank == 0)
            {
                printf("Attempting to generate %d structures in spacegroup number %d....\n",
                spg_num_structures, spg);
            }

        while( counter < spg_num_structures )
        {
            int verdict = 0; //for structure check
            long i = 0;       //counts attempts for spg
            //attempts for an spg.
            stop_flag = 0;
          
            for(; i < (long) max_attempts/total_ranks; i = i + BATCH_SIZE)

            {
                long j = 0;
                success_flag = 0;
                for(; j < BATCH_SIZE; j++)
                {
                    //generate
                    int result = generate_crystal(random_crystal,
                                                  mol,
                                                  volume,
                                                  Z,
                                                  Zp_max,
                                                  spg,
                                                  compatible_spg,
                                                  num_compatible_spg,
                                                  spg_index,
                                                  norm_dev,
                                                  angle_std);


                    //reset volume after volume attempts
                    if( (i+j) % vol_attempt == 0 && i+j != 0)
                    {
                        do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
                        if(my_rank == 0)
                            printf("#Rank %8d: Completed %ld attempts.\n", 0, (long)(i+j)*total_ranks);
                        //printf("fail count - %d / %d \n", fail_count, i);
                        fflush(stdout);
                    }

                    //alignment failure
                    if(!result)
                    {
                        //fail_count ++;
                        continue;
                    }

		    //check if molecules are too close with sr
                    verdict = check_structure_with_vdw_matrix(*random_crystal, vdw_matrix, dim1, dim2);
                    //if generation is successful
                    if (verdict == 1)
                    {
                        random_crystal->spg = spg;
                        break;
                    }
                }//end of GRAIN loop

                int found_poll[total_ranks];
                MPI_Gather(&verdict, 1, MPI_INT, &found_poll, 1, MPI_INT, 0, world_comm);
                if (my_rank == 0)
                {
                    //print the structure generated by root first to outfile
                    if(verdict)
                    {
                        if (counter < spg_num_structures)
                        {
                            print_crystal2file(random_crystal, out_file);
                            printf("#Rank %8d: Generation successful.\n", my_rank);
                            counter++;
                            success_flag = 1;
                            //int spglib_spg = detect_spg_using_spglib(random_crystal);
                            //printf("#SPGLIB detected space group = %d\n\n",
                            //                                   spglib_spg);
                            //printf("fail count = %d, total_attempts = %d \n", fail_count, i);
                            //fail_count = 0;

                        }
                        else
                            stop_flag = 1;

                    }
                    //get and print structures from other ranks
                    for(int rank = 1; rank < total_ranks; rank++)
                    {
                        if (found_poll[rank] == 1)
                        {
                            receive_xtal(world_comm, rank, random_crystal, total_atoms);
                            success_flag = 1;
                            if(counter < spg_num_structures)
                            {
                                print_crystal2file(random_crystal, out_file);
                                printf("#Rank %8d: Generation successful.\n", rank);
                                counter++;
                            }
                            else
                                stop_flag = 1;
                        }
                    }
                }
                //all other ranks send the crystal to root rank
                else
                {
                    if (verdict == 1)
                    {
                        send_xtal(world_comm, 0, random_crystal, total_atoms);
                        i = 0;
                    }
                }

                MPI_Bcast(&success_flag, 1, MPI_INT, 0, world_comm);
                MPI_Bcast(&stop_flag, 1, MPI_INT, 0, world_comm);

                // Atleast one structure was generated succesfully.
                // Reset counter and unit cell volume.
                if (success_flag)
                {
                    i = 0;
                    do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
                }

                // Enough structures were generated.
                if (stop_flag)
                    break;

            }//end of attempt loop

            //if max limit is reached if some rank hit the limit
            if (i >= max_attempts/total_ranks)
            {
                do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
                if (my_rank== 0)
                {
                    printf("**WARNING: generation failed for space group = %d "
                            "after %12ld attempts. \n",
                            spg,
                            max_attempts);
                    fflush(stdout);
                    //print_crystal(random_crystal);
                }
                counter = spg_num_structures + 1;
                MPI_Barrier(world_comm);
            }

            if(stop_flag)
                break;

        }//end of numof structures whileloop

        //move to next spacegroup
        counter = 0;
        spg_index++;
        if (my_rank == 0)
        {
             printf("#space group counter reset. Moving to next space group...\n\n");
             fflush(stdout);
        }
        MPI_Barrier(world_comm);

        //timing informatino
        time_t end_time = time(NULL);
        double elapsed = difftime (end_time, start_time);
        if (my_rank == 0)
        {
            printf("\nTIMING INFO:\n");
            printf("-----------------------------------------------------\n");
            print_time();
            printf("Time spent on space group %d: ~ %.0lf seconds \n", spg, elapsed);
            printf("-----------------------------------------------------\n\n");
        }

    }//end of spg while loop

    if(my_rank == 0)
        fclose(out_file);

    if(my_rank == 0)
    {
        print_time();
        printf("Generation completed.\nHave a nice day!!\n");
    }

}

            ////////////////////////////// Below for layer generation //////////////////////////////////////////

void mpi_generate_layer_with_vdw_cutoff_matrix(
	float *vdw_matrix,
	int dim1,
	int dim2,
	int num_structures,
	int Z,
	double volume_mean1,
	double volume_std1,
	double interface_area_mean1,
	double interface_area_std1,
	int volume_multiplier,
	double tol1,
	long max_attempts,
	char *spg_dist_type,
	float lattice_vector_2d_from_geo[2][3],
    int vol_attempt,
    int random_seed,
	MPI_Comm world_comm)
{
	float Zp_max=1;
    	float volume_mean = volume_mean1;
    	float volume_std = volume_std1;

	//added interface surface area mean and std
	float interface_area_mean = interface_area_mean1;
	float interface_area_std = interface_area_std1;
    	TOL = tol1;
	int SET_INTERFACE_AREA = 0;

	// if dont want to manually set it, set interface_area_mean = 0 in ui.conf
	if (interface_area_mean != 0)
	{
		SET_INTERFACE_AREA = 1;
	}

    	if (dim1 != dim2)
	{printf("***ERROR:vdw cutoff matrix is not square***\n"); exit(0);}

	//Initialise MPI
	int total_ranks;
    	MPI_Comm_size(world_comm, &total_ranks);
    	int my_rank;
   	MPI_Comm_rank(world_comm, &my_rank);

	//random number seeding
	srand((unsigned int)time(NULL));
	int seed_shift = rand()% 1061 + 7 ;
	int stop_flag = 0;	// to stop threads if limit is reached
	int success_flag = 0;
	int counter = 0;	//counts number of structures
	int spg_index = 0;	//space group to be generated
	FILE *out_file;		//file to output geometries

	/*
	global_lattice_vector_2d[0][0]=lattice_vector_2d_from_geo[0][0];
	global_lattice_vector_2d[0][1]=lattice_vector_2d_from_geo[0][1];
	global_lattice_vector_2d[0][2]=lattice_vector_2d_from_geo[0][2];
	global_lattice_vector_2d[1][0]=lattice_vector_2d_from_geo[1][0];
	global_lattice_vector_2d[1][1]=lattice_vector_2d_from_geo[1][1];
	global_lattice_vector_2d[1][2]=lattice_vector_2d_from_geo[1][2];
 	*/
	if (my_rank == 0)
	{
		out_file = fopen("geometry.out","w");
		if(!out_file)		//check permissions
		{
			printf("***ERROR: cannot create geometry.out \n");
			exit(0);
		}
	}
    else
    {
        out_file = NULL;
    }

	//random number seeding, different seeds for different threads
	seed = (unsigned int*)malloc(sizeof(unsigned int)); //seed for uniform gen
	seed2 = (unsigned int*)malloc(sizeof(unsigned int)); //seed for random
    	//using random seed from user
	//*seed = (unsigned int)abs(my_rank*7 + random_seed);  //some random seed private for each threads
    	//*seed2 = (unsigned int)abs(my_rank*17 + random_seed);
    	//does not use random seed from user
    	*seed += my_rank*7 + seed_shift*13; //some random seed private for each threads
	*seed2 = my_rank*17 + seed_shift*11;
	init_genrand(abs(*seed));

	//storing information for compatible space groups
	COMPATIBLE_SPG compatible_spg[230];
	int num_compatible_spg = 0;

	//variable declarartion
	molecule *mol = (molecule*)malloc(sizeof(molecule));//store molecule
	crystal *random_crystal = (crystal*)malloc(sizeof(crystal));//dummy crystal
	//float volume_std;	//standard dev for volumes
	//float volume_mean;	//mean volume
	float sr = -1;			//specific radius proportion for structure checks
						//see paper for definition
	//float Zp_max;		//Z'' . not implemented
	float volume;		//random volume used of generation
	//int Z;				//multiplicity of general position
	//int num_structures;	//num of structures
	int spg;			//space group attempted
	//long max_attempts;	//max attempts per space group

	//read input from file, read molecular structure from geometry,Z, Zp
	if(my_rank == 0)
	{
		int len = 100;
		char name[len];
		gethostname(name, len);
		printf("PARALLELIZATION INFO:\n");
		printf("---------------------------\n");
		printf("Using MPI for parallelization.\n");
		printf("cgenarris is running on %d processes. \n", total_ranks);
		printf("Master Rank host name: %s \n", name);
		print_time();
	}

	if(my_rank == 0)
	{
		printf("---------------------------\n");
		printf("\n");
	}
	MPI_Barrier(world_comm);

    read_geometry(mol, "geometry.in");				//read molecule from geometry.in

	//recenter molecule to origin
	recenter_molecule(mol);

	if (my_rank == 0)
	{
		print_input_geometry(mol);
		print_input_settings_layer(&num_structures,
							 &Z,
							 &Zp_max,
							 &volume_mean,
							 &volume_std,
                             				 &interface_area_mean,
                             				 &interface_area_std,
                             				 &volume_multiplier,
							 &sr,
                             				 lattice_vector_2d_from_geo,
							 &max_attempts,
							 spg_dist_type,
                             				 &vol_attempt,
                             				 &random_seed);
	}

	MPI_Barrier(world_comm);

    	//inititalise volume
    	do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
	int N = mol->num_of_atoms;
	allocate_xtal(random_crystal, ZMAX, N); //allcate memory
	random_crystal->num_atoms_in_molecule = mol->num_of_atoms;
	int total_atoms = mol->num_of_atoms * Z;

	//initialize an array to store substrate combos
	int max_num_combo = 100000000;
	int *all_substrate_combo= (int *) malloc(max_num_combo * 4 * sizeof (int));   //max 100000000 combinations
	int num_combo = generate_substrate_lattice_combs(all_substrate_combo,lattice_vector_2d_from_geo,volume,
							 MAX_ANGLE,MIN_ANGLE);



	if(my_rank == 0)
	{
		printf("COMPATIBLE LAYER GROUP INFO:\n");
		printf("-----------------------------------------------------\n");
		printf("Detecting compatible layer groups"
			   " from molecule symmetries...\n");
	}

	//compatible space groups are stored in compatible_spg. see the
	//object definition in spg_generation.h
	//every thread has its own copy
	find_compatible_lg_positions(mol,
				     Z,
				     compatible_spg,
				     &num_compatible_spg,
				     lattice_vector_2d_from_geo,
				     volume,
				     my_rank+1);

	if(my_rank == 0)
	{
		printf("\n");
		printf("Total compatible layer groups = %d\n", num_compatible_spg);
		//printf("Number of molecular axes = %d \n", num_axes);
		printf("-----------------------------------------------------\n\n");
		printf("Starting generation...\n");
		print_time();
		printf("\n");
		sleep(1);
	}

	MPI_Barrier(world_comm); // wait for other friends to join

	/*deprecated
	//find allowed space groups for general position (deprecated)
	//int allowed_spg[230];
	//int num_allowed_spg = 0;
	//find_allowed_spg(allowed_spg, &num_allowed_spg, Z);
	//printf("allowed %d \n", num_compatible_spg);
	*/


	while( spg_index < num_compatible_spg )
	{
		MPI_Barrier(world_comm);
		//counter counts the number of structures generated
		//if (spg_index > 44)
		//	break;
		//spg_index = 44;
		spg = compatible_spg[spg_index].spg; //pick a spg

		//time information
		time_t start_time = time(NULL);
		//fail_count = 0;

		// Get number of structures to be generated in this spg
		// skip to next spg if this is zero
		int spg_num_structures = find_num_structure_for_lg(num_structures, spg_dist_type, spg, Z);
		if( !spg_num_structures )
		{
			spg_index++;
			continue;
		}
		//print attempted space group
		if (my_rank == 0)
			{printf("Attempting to generate %d structures in layergroup number %d....\n", spg_num_structures, spg);}

		while( counter < spg_num_structures )
		{
			int verdict = 0; //for structure check
			int i = 0; 		 //counts attempts for spg
			//attempts for an spg.
			stop_flag = 0;
			for(; i < max_attempts/total_ranks; i = i + BATCH_SIZE)
			{
				int j = 0;
				success_flag = 0;
				for(; j < BATCH_SIZE; j++)
				{
					//generate

					int result = generate_layer_crystal(random_crystal,
									 mol,
									 volume,
									 Z,
									 Zp_max,
									 spg,
									 compatible_spg,
									 num_compatible_spg,
									 spg_index,
									 lattice_vector_2d_from_geo,
									 all_substrate_combo,
									 num_combo,
									 interface_area_mean,
									 interface_area_std,
									 volume_multiplier,
									 SET_INTERFACE_AREA);

					//printf("I am after generate_crystal\n");

					//reset volume after volume attempts
					if( (i+j) % vol_attempt == 0 && i+j != 0)
					{
						do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
						if(my_rank == 0)
							printf("#Rank %8d: Completed %d attempts.\n", 0, (i+j)*total_ranks);
						//printf("fail count - %d / %d \n", fail_count, i);
						//print_crystal(random_crystal);
						//fflush(stdout);
					}

					if (result ==2)
					{
						//printf("I am in result==2");
						//fflush(stdout);
						break;
					}

					//alignment failure
					if(!result)
					{
						//fail_count ++;
						continue;
					}
					//check if molecules are too close with sr
					verdict = check_structure_with_vdw_matrix(*random_crystal, vdw_matrix, dim1, dim2);
					//printf("I am after check_structure_with_vdw_matrix\n");
					//fflush(stdout);
										//if generation is successful
					if (verdict == 1 )
					{
						random_crystal->spg = spg;
						break;
					}


				}//end of GRAIN loop

				int found_poll[total_ranks];
				//printf("verdict = %d\n", verdict);
				MPI_Gather(&verdict, 1, MPI_INT, &found_poll, 1, MPI_INT, 0, world_comm);
				if (my_rank == 0)
				{
					//print the structure generated by root first to outfile
					if(verdict)
					{
						if (counter < spg_num_structures)
						{
                            //add a layer group check here
                            //check_layer_group(random_crystal);
							print_layer2file(random_crystal, out_file);
							printf("#Rank %8d: Generation successful.\n", my_rank);
							counter++;
							success_flag = 1;
							//int spglib_spg = detect_spg_using_spglib(random_crystal);
							//printf("#SPGLIB detected space group = %d\n\n",
							//				                     spglib_spg);
							//printf("fail count = %d, total_attempts = %d \n", fail_count, i);
							//fail_count = 0;

						}
						else
							stop_flag = 1;

					}
					//get and print structures from other ranks
					for(int rank = 1; rank < total_ranks; rank++)
					{
						if (found_poll[rank] == 1)
						{
							receive_xtal(world_comm, rank, random_crystal, total_atoms);
							//print_crystal(random_crystal);
							success_flag = 1;
							if(counter < spg_num_structures)
							{
								print_layer2file(random_crystal, out_file);
								printf("#Rank %8d: Generation successful.\n", rank);
								counter++;
							}
							else
								stop_flag = 1;
						}
					}
				}
				//all other ranks send the crystal to root rank
				else
				{
					if (verdict == 1)
					{
						send_xtal(world_comm, 0, random_crystal, total_atoms);
						i = 0;
					}
				}

				MPI_Bcast(&success_flag, 1, MPI_INT, 0, world_comm);
				MPI_Bcast(&stop_flag, 1, MPI_INT, 0, world_comm);

				if (success_flag)
				{
					i = 0;
					do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
				}
				if (stop_flag)
					break;


			}//end of attempt loop

			//if max limit is reached if some rank hit the limit
			if (i >= max_attempts/total_ranks)
			{
				do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
				if (my_rank== 0)
				{
					printf("**WARNING: generation failed for layer group = %d "
							"after %12ld attempts. \n",
							spg,
							max_attempts);
					fflush(stdout);
					//print_crystal(random_crystal);
				}
				counter = spg_num_structures + 1;
				MPI_Barrier(world_comm);
			}

			if(stop_flag)
				break;

		}//end of numof structures whileloop

		//move to next spacegroup
		counter = 0;
		spg_index++;
		if (my_rank == 0)
		{
			 printf("#layer group counter reset. Moving to next layer group...\n\n");
			 fflush(stdout);
		}
		MPI_Barrier(world_comm);

		//timing informatino
		time_t end_time = time(NULL);
		double elapsed = difftime (end_time, start_time);
		if (my_rank == 0)
		{
			printf("\nTIMING INFO:\n");
			printf("-----------------------------------------------------\n");
			print_time();
			printf("Time spent on layer group %d: ~ %.0lf seconds \n", spg, elapsed);
			printf("-----------------------------------------------------\n\n");
		}

	}//end of spg while loop

	if(my_rank == 0)
		fclose(out_file);

	if(my_rank == 0)
	{
		print_time();
		printf("Generation completed.\nHave a nice day!!\n");
	}

	free(all_substrate_combo);

}



void send_xtal(MPI_Comm comm, int destination, crystal* xtal, int total_atoms)
{
    //2d array memory need not be continous. copy to 1d then send.
    float temp[9] = {xtal->lattice_vectors[0][0], xtal->lattice_vectors[0][1],
                     xtal->lattice_vectors[0][2], xtal->lattice_vectors[1][0],
                     xtal->lattice_vectors[1][1], xtal->lattice_vectors[1][2],
                     xtal->lattice_vectors[2][0], xtal->lattice_vectors[2][1],
                     xtal->lattice_vectors[2][2] };
    MPI_Send(temp, 9, MPI_FLOAT , destination, 1, comm);
    MPI_Send(xtal->Xcord, total_atoms, MPI_FLOAT , destination, 2, comm);
    MPI_Send(xtal->Ycord, total_atoms, MPI_FLOAT , destination, 3, comm);
    MPI_Send(xtal->Zcord, total_atoms, MPI_FLOAT , destination, 4, comm);
    MPI_Send(xtal->atoms, 2*total_atoms, MPI_CHAR , destination, 5, comm);
    MPI_Send(&(xtal->spg), 1, MPI_INT , destination, 6, comm);
    MPI_Send(&(xtal->wyckoff_position), 1, MPI_INT , destination, 7, comm);
    MPI_Send(&(xtal->num_atoms_in_molecule), 1, MPI_INT , destination, 8, comm);
    MPI_Send(&(xtal->Z), 1, MPI_INT , destination, 9, comm);
    MPI_Send(&(xtal->Zp), 1, MPI_INT , destination, 10, comm);
    MPI_Send(xtal->com_positions, 3, MPI_FLOAT , destination, 11, comm);
    MPI_Send(xtal->euler_angles, 3, MPI_FLOAT , destination, 12, comm);
}

void receive_xtal(MPI_Comm comm, int source, crystal* xtal, int total_atoms)
{
    MPI_Status status;
    float temp[9];
    MPI_Recv(temp, 9, MPI_FLOAT, source, 1, comm, &status);
    MPI_Recv(xtal->Xcord, total_atoms, MPI_FLOAT, source, 2, comm, &status);
    MPI_Recv(xtal->Ycord, total_atoms, MPI_FLOAT, source, 3, comm, &status);
    MPI_Recv(xtal->Zcord, total_atoms, MPI_FLOAT, source, 4, comm, &status);
    MPI_Recv(xtal->atoms, 2*total_atoms, MPI_CHAR, source, 5, comm, &status);
    MPI_Recv(&(xtal->spg), 1, MPI_INT, source, 6, comm, &status);
    MPI_Recv(&(xtal->wyckoff_position), 1, MPI_INT, source, 7, comm, &status);
    MPI_Recv(&(xtal->num_atoms_in_molecule), 1, MPI_INT, source, 8, comm, &status);
    MPI_Recv(&(xtal->Z), 1, MPI_INT, source, 9, comm, &status);
    MPI_Recv(&(xtal->Zp), 1, MPI_INT, source, 10, comm, &status);

    MPI_Recv(xtal->com_positions, 3, MPI_FLOAT, source, 11, comm, &status);
    MPI_Recv(xtal->euler_angles, 3, MPI_FLOAT, source, 12, comm, &status);
    
    xtal->lattice_vectors[0][0] = temp[0];
    xtal->lattice_vectors[0][1] = temp[1];
    xtal->lattice_vectors[0][2] = temp[2];

    xtal->lattice_vectors[1][0] = temp[3];
    xtal->lattice_vectors[1][1] = temp[4];
    xtal->lattice_vectors[1][2] = temp[5];

    xtal->lattice_vectors[2][0] = temp[6];
    xtal->lattice_vectors[2][1] = temp[7];
    xtal->lattice_vectors[2][2] = temp[8];
}

int num_compatible_spacegroups(int Z, double tolerance)
{
    //set global variable tolerance
    TOL = tolerance;

    COMPATIBLE_SPG compatible_spg[230];
    int num_compatible_spg = 0;
    int thread_num = 1;
    molecule *mol = (molecule*)malloc(sizeof(molecule));

    //read geometry from geometry.in
    read_geometry(mol, "geometry.in");

    find_compatible_spg_positions(mol,
                                  Z,
                                  compatible_spg,
                                  &num_compatible_spg,
                                  thread_num);

    free(mol);

    return num_compatible_spg;

}

///////////////// for layer group /////////////////////

int num_compatible_layergroups(int Z, double tolerance,float volume,float lattice_vector_2d_from_geo[2][3])
{
	//set global variable tolerance
	TOL = tolerance;

	COMPATIBLE_SPG compatible_spg[230];
	int num_compatible_lg = 0;
	int thread_num = 1;
	molecule *mol = (molecule*)malloc(sizeof(molecule));

	//read geometry from geometry.in
	read_geometry(mol, "geometry.in");
    //float global_lattice_vector_2d[2][3] = {{5.7675159999999996,0.0000000000000000,0.0000000000000000},
											//{2.8837579999999998,4.9948160000000001,0.0000000000000000}};
	//hard code global lattice vector here for now, delete latter
	find_compatible_lg_positions(mol,
				     Z,
				     compatible_spg,
				     &num_compatible_lg,
			             lattice_vector_2d_from_geo,
				     volume,
			             thread_num);

	free(mol);
	return num_compatible_lg;

}

void create_crystal_from_array(crystal *xtal, double lattice_vector[3][3],
    double *Xc,int total_atoms1, double *Yc,int total_atoms2,
    double *Zc, int total_atoms3, char* atoms, int total_atoms,
    int Z, int spg)
{

    xtal->spg = spg;
    xtal->Z = Z;

    if (total_atoms != total_atoms1 || total_atoms != total_atoms2 || total_atoms != total_atoms3)
    printf("***ERR0R:sizes of arrays doesnot match");

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

}

void create_array_from_crystal(crystal *xtal, double lattice_vector[3][3],
    double *Xc,int total_atoms1, double *Yc,int total_atoms2,
    double *Zc, int total_atoms3, char* atoms, int total_atoms,
    int Z, int spg)
{
    spg = xtal->spg;
    Z = xtal->Z;

    if (total_atoms != total_atoms1 || total_atoms != total_atoms2 || total_atoms != total_atoms3)
    printf("***ERR0R:sizes of arrays doesnot match");


    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
	lattice_vector[i][j] = xtal->lattice_vectors[i][j];
    }
    for(int i = 0; i < total_atoms; i++)
    {
        Xc[i] = xtal->Xcord[i];
        Yc[i] = xtal->Ycord[i];
        Zc[i] = xtal->Zcord[i];
        atoms[2*i]   = xtal->atoms[2*i];
        atoms[2*i+1] = xtal->atoms[2*i+1];
    }
}

