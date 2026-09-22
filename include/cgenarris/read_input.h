#ifndef _READ_INPUT_H_
#define _READ_INPUT_H_

#include <stdio.h>
#include <stdlib.h>

#include "cgenarris/molecule.h"
#include "cgenarris/input_settings.h"

/*
Fills set with the defaults below, then reads "key value" lines from
path. Blank lines and lines starting with '#' are ignored; keys are
case sensitive; an unknown key, a malformed value or a missing required
key is an error.

Keys (unit) [default]:
    generation_type         crystal | layer | asu           [crystal]
    number_of_structures    per space group, or total ASUs  (required)
    max_attempts            per space group, or per ASU     [100000]
    random_seed             0 = time based                  [0]
    tolerance               special-position tolerance      [0.1]
  crystal / layer
    Z                       molecules per cell              (required)
    volume_mean             Angstrom^3                      (required)
    volume_std              Angstrom^3                      (required)
    volume_attempts         redraw volume after N failures  [100000]
    sr                      specific radius                 [0.85]
    lattice_norm_dev                                        [0.4]
    lattice_angle_std       degrees                         [8]
    spg_distribution_type   standard|uniform|csd|chiral|racemic [standard]
    rigid_press             0 | 1                           [0]
  layer only
    interface_area_mean     Angstrom^2                      [0]
    interface_area_std      Angstrom^2                      [0]
    volume_multiplier                                       [3]
    lattice_vector_a        x y z, Angstrom                 [0 0 0]
    lattice_vector_b        x y z, Angstrom                 [0 0 0]
  multi-component / asu (molecules are read from geometry_<i>.in)
    molecule_types          number of distinct molecules    [1]
    stoichiometry           one integer per molecule type   [1 ...]
    asu_sr_min              sr window lower bound           [0.75]
    asu_sr_max              sr window upper bound           [1.3]
    asu_output_file         output path, geometry.out format [asu.out]

Returns:
    0 on success, -1 on error (message printed to stderr).
*/
int read_settings(Settings *set, const char *path);

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
