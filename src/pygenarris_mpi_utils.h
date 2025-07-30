#ifndef _PYGENARRIS_MPI_UTILS_H_
#define _PYGENARRIS_MPI_UTILS_H_

#include "mpi.h"

#include "input_settings.h"
#include "molecule.h"

#define NO_STOP        0
#define ENOUGH_STOP    1
#define ATTEMPTS_STOP  2


void print_time(void);
FILE* open_output_file(int my_rank);
void init_random_seed(unsigned int *seed, unsigned int *seed2,
                      int random_seed, int rank);
void recenter_molecules(molecule* mol, int mol_types);
float draw_volume(float volume_mean, float volume_std);
void get_n_atoms_in_mol(int *n_atoms_in_mol, molecule *mol, int n_mol_types);
void print_allowed_spg(int *allowed_spg, int num_spg);
int check_stop_condition(int struct_counter, int max_struct,
                         long attempt, long max_attempt);
void print_spg_end(double elapsed, int struct_counter, int spg);
void print_exit();


#endif //  pygenarris_mpi_utils.h
