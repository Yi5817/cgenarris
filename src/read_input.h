#ifndef _READ_INPUT_H_
#define _READ_INPUT_H_

#include <stdio.h>
#include <stdlib.h>

#include "molecule.h"
#include "input_settings.h"

void read_control(int* num_structures, int* Z, float* Zp_max,
                 float* volume_mean,float* volume_std,
                 float *sr,long *max_attempts, char *spg_dist_type,
                 int *vol_attempts,int *random_seed,
                 char *crystal_generation,
                 float* interface_area_mean,float* interface_area_std,
                 int* volume_multiplier,
                 float lattice_vector_2d[2][3], float* norm_dev,
		 float* angle_std, int **stoic, int *mol_types,
		 int *rigid_press);

void read_geometry(molecule* mol, char* filename);
void read_molecules(molecule *mol, int mol_types);

void print_input_geometry(molecule* mol);
void print_input_geometries(molecule *mol, int mol_types);
void print_molecule(molecule *mol);
void print_input_settings(Settings set);

void print_input_settings_layer(int* num_structures,
                          int* Z,
                          float* Zp_max,
                          float* volume_mean,
                          float* volume_std,
                          float* interface_area_mean,
                          float* interface_area_std,
                          int* volume_multiplier,
                          float *sr,
                          float global_lattice_vector_2d[2][3],
                          long *max_attempts,
                          char * spg_dist_type,
                          int *vol_attempt,
                          int *random_seed);

#endif
