#ifndef _INPUT_SETTINGS_H_
#define _INPUT_SETTINGS_H_

//generation types
#define COCRYSTAL 1
#define CRYSTAL   2
#define LAYER     3

typedef struct
{
    int num_structures;
    int Z;
    int vol_attempts;
    int random_seed;
    int *stoic;
    int n_mol_types;
    int generation_type;
    int rigid_press;
    long max_attempts;
    char *spg_dist_type;
    float vol_mean;
    float vol_std;
    float norm_dev;
    float angle_std;
    float sr;
    float *vdw_matrix;

}Settings;


#endif
