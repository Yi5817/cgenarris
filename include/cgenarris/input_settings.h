#ifndef _INPUT_SETTINGS_H_
#define _INPUT_SETTINGS_H_

//generation types
#define CRYSTAL 1
#define LAYER   2
#define ASU     3

// Largest number of distinct molecule types in a multi-component run.
#define MAX_MOL_TYPES 8

/*
All user-facing settings of a cgenarris run. Filled from control.in by
read_settings() (see read_input.h for keys, units and defaults), or field
by field from the Python bindings.
*/
typedef struct
{
    // Common
    int generation_type;          // CRYSTAL, LAYER or ASU
    int num_structures;           // crystal/layer: per space group; asu: total
    long max_attempts;            // crystal/layer: per space group; asu: per unit
    int random_seed;              // 0 = time based
    float tol;                    // numerical tolerance for special positions

    // Crystal and layer
    int Z;                        // molecules (4.0: asymmetric units) per cell
    float vol_mean;               // Angstrom^3
    float vol_std;                // Angstrom^3
    int vol_attempts;             // failed attempts before the volume is redrawn
    float sr;                     // specific radius for the vdW check
    float norm_dev;               // lattice principal-component deviation
    float angle_std;              // lattice angle std, degrees
    char spg_dist_type[16];       // standard | uniform | csd | chiral | racemic
    int rigid_press;
    float *vdw_matrix;            // set by the generator, not parsed

    // Layer only
    float interface_area_mean;    // Angstrom^2
    float interface_area_std;     // Angstrom^2
    int volume_multiplier;
    float lattice_vector_2d[2][3];

    // Multi-component / asymmetric unit
    int n_mol_types;              // 0 = single molecule from geometry.in
    int stoic[MAX_MOL_TYPES];     // copies of each molecule type per unit
    float asu_sr_min;             // sr window lower bound (exclusive)
    float asu_sr_max;             // sr window upper bound (exclusive)
    char asu_output_file[256];    // output path, geometry.out block format

}Settings;


#endif
