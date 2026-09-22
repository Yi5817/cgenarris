#include "cgenarris/read_input.h"
#include "molecule_utils.h"
#include "cgenarris/input_settings.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


float TOL;

// ---- control.in parser -----------------------------------------------------

#define LINE_MAX_LEN 1024

static void settings_defaults(Settings *set)
{
    memset(set, 0, sizeof(*set));
    set->generation_type = CRYSTAL;
    set->num_structures = -1;          // required
    set->max_attempts = 100000;
    set->random_seed = 0;
    set->tol = 0.1;
    set->Z = -1;                       // required for crystal/layer
    set->vol_mean = -1;                // required for crystal/layer
    set->vol_std = -1;                 // required for crystal/layer
    set->vol_attempts = 100000;
    set->sr = 0.85;
    set->norm_dev = 0.4;
    set->angle_std = 8;
    strcpy(set->spg_dist_type, "standard");
    set->rigid_press = 0;
    set->vdw_matrix = NULL;
    set->interface_area_mean = 0;
    set->interface_area_std = 0;
    set->volume_multiplier = 3;
    set->n_mol_types = 1;
    for(int i = 0; i < MAX_MOL_TYPES; i++)
        set->stoic[i] = 0;             // 0 = not given; filled with 1s later
    set->asu_sr_min = 0.75;
    set->asu_sr_max = 1.3;
    strcpy(set->asu_output_file, "asu.out");
}

// Parses one integer token; returns 0 on success.
static int parse_int(const char *key, const char *tok, int *out)
{
    char *end;
    long v = strtol(tok, &end, 10);
    if(tok[0] == '\0' || *end != '\0')
    {
        fprintf(stderr, "***ERROR: control.in: %s expects an integer, got '%s'\n", key, tok);
        return -1;
    }
    *out = (int)v;
    return 0;
}

static int parse_long(const char *key, const char *tok, long *out)
{
    char *end;
    long v = strtol(tok, &end, 10);
    if(tok[0] == '\0' || *end != '\0')
    {
        fprintf(stderr, "***ERROR: control.in: %s expects an integer, got '%s'\n", key, tok);
        return -1;
    }
    *out = v;
    return 0;
}

static int parse_float(const char *key, const char *tok, float *out)
{
    char *end;
    double v = strtod(tok, &end);
    if(tok[0] == '\0' || *end != '\0')
    {
        fprintf(stderr, "***ERROR: control.in: %s expects a number, got '%s'\n", key, tok);
        return -1;
    }
    *out = (float)v;
    return 0;
}

/*
Applies one "key value..." line to set. `value` is the first value token;
further tokens are fetched with strtok(NULL, ...). Returns 0 on success.
*/
static int apply_setting(Settings *set, const char *key, char *value)
{
    if(!value)
    {
        fprintf(stderr, "***ERROR: control.in: %s has no value\n", key);
        return -1;
    }

    if(!strcmp(key, "generation_type"))
    {
        if(!strcmp(value, "crystal"))    set->generation_type = CRYSTAL;
        else if(!strcmp(value, "layer")) set->generation_type = LAYER;
        else if(!strcmp(value, "asu"))   set->generation_type = ASU;
        else
        {
            fprintf(stderr, "***ERROR: control.in: generation_type must be "
                    "crystal, layer or asu, got '%s'\n", value);
            return -1;
        }
        return 0;
    }
    if(!strcmp(key, "number_of_structures")) return parse_int(key, value, &set->num_structures);
    if(!strcmp(key, "max_attempts"))         return parse_long(key, value, &set->max_attempts);
    if(!strcmp(key, "random_seed"))          return parse_int(key, value, &set->random_seed);
    if(!strcmp(key, "tolerance"))            return parse_float(key, value, &set->tol);

    if(!strcmp(key, "Z"))                    return parse_int(key, value, &set->Z);
    if(!strcmp(key, "volume_mean"))          return parse_float(key, value, &set->vol_mean);
    if(!strcmp(key, "volume_std"))           return parse_float(key, value, &set->vol_std);
    if(!strcmp(key, "volume_attempts"))
    {
        int v;
        if(parse_int(key, value, &v))
            return -1;
        if(v != 0)                           // 0 keeps the default
            set->vol_attempts = v;
        return 0;
    }
    if(!strcmp(key, "sr"))                   return parse_float(key, value, &set->sr);
    if(!strcmp(key, "lattice_norm_dev") || !strcmp(key, "norm_dev"))
        return parse_float(key, value, &set->norm_dev);
    if(!strcmp(key, "lattice_angle_std") || !strcmp(key, "angle_std"))
        return parse_float(key, value, &set->angle_std);
    if(!strcmp(key, "rigid_press"))          return parse_int(key, value, &set->rigid_press);
    if(!strcmp(key, "spg_distribution_type"))
    {
        if(strcmp(value, "standard") && strcmp(value, "uniform") &&
           strcmp(value, "chiral") && strcmp(value, "racemic") && strcmp(value, "csd"))
        {
            fprintf(stderr, "***ERROR: control.in: bad value of "
                    "spg_distribution_type '%s'\n", value);
            return -1;
        }
        strncpy(set->spg_dist_type, value, sizeof(set->spg_dist_type) - 1);
        return 0;
    }

    if(!strcmp(key, "interface_area_mean"))  return parse_float(key, value, &set->interface_area_mean);
    if(!strcmp(key, "interface_area_std"))   return parse_float(key, value, &set->interface_area_std);
    if(!strcmp(key, "volume_multiplier"))    return parse_int(key, value, &set->volume_multiplier);
    if(!strcmp(key, "lattice_vector_a") || !strcmp(key, "lattice_vector_b"))
    {
        float *vec = set->lattice_vector_2d[key[strlen(key) - 1] - 'a'];
        if(parse_float(key, value, vec))
            return -1;
        for(int i = 1; i < 3; i++)
        {
            char *tok = strtok(NULL, " \t\r\n");
            if(!tok || parse_float(key, tok, vec + i))
            {
                fprintf(stderr, "***ERROR: control.in: %s expects three numbers\n", key);
                return -1;
            }
        }
        return 0;
    }

    if(!strcmp(key, "molecule_types"))
    {
        if(parse_int(key, value, &set->n_mol_types))
            return -1;
        if(set->n_mol_types < 1 || set->n_mol_types > MAX_MOL_TYPES)
        {
            fprintf(stderr, "***ERROR: control.in: molecule_types must be in "
                    "[1, %d], got %d\n", MAX_MOL_TYPES, set->n_mol_types);
            return -1;
        }
        return 0;
    }
    if(!strcmp(key, "stoichiometry") || !strcmp(key, "stochiometry"))
    {
        int n = 0;
        char *tok = value;
        while(tok)
        {
            if(n == MAX_MOL_TYPES)
            {
                fprintf(stderr, "***ERROR: control.in: stoichiometry has more "
                        "than %d entries\n", MAX_MOL_TYPES);
                return -1;
            }
            if(parse_int(key, tok, set->stoic + n))
                return -1;
            if(set->stoic[n] < 1)
            {
                fprintf(stderr, "***ERROR: control.in: stoichiometry entries "
                        "must be >= 1, got %d\n", set->stoic[n]);
                return -1;
            }
            n++;
            tok = strtok(NULL, " \t\r\n");
        }
        for(int i = n; i < MAX_MOL_TYPES; i++)
            set->stoic[i] = 0;         // marks the end of the given entries
        return 0;
    }
    if(!strcmp(key, "asu_sr_min"))           return parse_float(key, value, &set->asu_sr_min);
    if(!strcmp(key, "asu_sr_max"))           return parse_float(key, value, &set->asu_sr_max);
    if(!strcmp(key, "asu_output_file"))
    {
        if(strlen(value) >= sizeof(set->asu_output_file))
        {
            fprintf(stderr, "***ERROR: control.in: asu_output_file is too long\n");
            return -1;
        }
        strcpy(set->asu_output_file, value);
        return 0;
    }

    fprintf(stderr, "***ERROR: control.in: unknown key '%s'\n", key);
    return -1;
}

// Cross-field checks after every line has been read.
static int validate_settings(Settings *set)
{
    const char *problem = NULL;

    if(set->num_structures < 0)
        problem = "number_of_structures is required";
    else if(set->max_attempts < 1)
        problem = "max_attempts must be >= 1";
    else if(set->generation_type == ASU)
    {
        int n_given = 0;
        while(n_given < MAX_MOL_TYPES && set->stoic[n_given] > 0)
            n_given++;
        if(n_given == 0)                       // key omitted: one copy of each
            for(int i = 0; i < set->n_mol_types; i++)
                set->stoic[i] = 1;
        else if(n_given != set->n_mol_types)
            problem = "stoichiometry must list one entry per molecule type";

        if(!problem && !(set->asu_sr_min >= 0 && set->asu_sr_min < set->asu_sr_max))
            problem = "asu_sr_min must satisfy 0 <= asu_sr_min < asu_sr_max";
    }
    else
    {
        if(set->Z < 1)
            problem = "Z is required and must be >= 1";
        else if(set->vol_mean <= 0)
            problem = "volume_mean is required and must be > 0";
        else if(set->vol_std < 0)
            problem = "volume_std is required and must be >= 0";
        else if(set->sr <= 0)
            problem = "sr must be > 0";
    }

    if(problem)
    {
        fprintf(stderr, "***ERROR: control.in: %s\n", problem);
        return -1;
    }
    return 0;
}

int read_settings(Settings *set, const char *path)
{
    settings_defaults(set);

    FILE *fileptr = fopen(path, "r");
    if(!fileptr)
    {
        fprintf(stderr, "***ERROR: cannot open %s\n", path);
        return -1;
    }

    char line[LINE_MAX_LEN];
    int line_no = 0;
    int status = 0;
    while(status == 0 && fgets(line, sizeof(line), fileptr))
    {
        line_no++;
        char *key = strtok(line, " \t\r\n");
        if(!key || key[0] == '#')
            continue;
        char *value = strtok(NULL, " \t\r\n");
        if(apply_setting(set, key, value))
        {
            fprintf(stderr, "***ERROR: %s line %d could not be read\n", path, line_no);
            status = -1;
        }
    }
    fclose(fileptr);

    if(status == 0)
        status = validate_settings(set);
    if(status == 0)
        TOL = set->tol;
    return status;
}


void read_geometry(molecule* mol, char* filename)
{
    FILE *fileptr;
    size_t len = 0;
    int read;
    char *line = NULL;
    char *sub_line = NULL;
    int i = 0;
    int atom_count = 0;

    //find_number_of atoms
    fileptr = fopen(filename,"r");
    //check if file exits
    if(!fileptr)
    {
        printf("***ERROR: no %s file \n", filename);
        exit(EXIT_FAILURE);
    }

     while ((read = getline(&line, &len, fileptr)) != -1)
    {
        if (strstr(line, "#") != NULL)
            continue;

        sub_line=strtok(line," ");
        //printf("%s \n" , sub_line);
        if(strcmp(sub_line, "atom") == 0)
            atom_count++;
        else
            continue;
    }

    //printf("Total number of atoms in molecule = %d\n", atom_count);
    int N = atom_count;
    //memory allocation
    (*mol).atoms = (char *)malloc(2*N*sizeof(char));
    (*mol).X = (float *)malloc(N*sizeof(float));
    (*mol).Y = (float *)malloc(N*sizeof(float));
    (*mol).Z = (float *)malloc(N*sizeof(float));

    fseek(fileptr, 0, SEEK_SET);
    while ((read = getline(&line, &len, fileptr)) != -1)
    {
        //printf("%s", line);
        if (strstr(line, "#") != NULL)
            continue;

        sub_line=strtok(line," ");
        if(strcmp(sub_line, "atom") == 0)
            atom_count++;
        else
            continue;
        sub_line=strtok(NULL," ");
        (*mol).X[i]=atof(sub_line);
        sub_line=strtok(NULL," ");
        (*mol).Y[i]=atof(sub_line);
        sub_line=strtok(NULL," ");
        (*mol).Z[i]=atof(sub_line);
        sub_line=strtok(NULL," ");
        (*mol).atoms[2*i]=*sub_line;
        if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
            (*mol).atoms[2*i+1]=' ';
        else
            (*mol).atoms[2*i+1]=*(sub_line+1);
        i++;

    }
    fclose(fileptr);
   // printf("atoms = %d \n", atom_count);
    mol->num_of_atoms = N;

}

void read_molecules(molecule *mol, int mol_types)
{
    for(int i = 0; i < mol_types; i++)
    {
        char filename[25];
        sprintf(filename, "geometry_%d.in", i);
        read_geometry(mol + i, filename);
    }

}

void print_input_geometries(molecule *mol, int mol_types)
{
    printf("Total number of molecule geometries = %d\n\n", mol_types);
    for(int i = 0; i < mol_types; i++)
    {
        printf("MOLECULE GEOMETRY %d:\n", i);
        printf("-----------------------------\n");
        print_molecule(mol + i);
        printf("-----------------------------\n\n");
    }
}


void print_input_settings(Settings set)
{
    printf("INPUT SETTINGS:\n");
    printf("-----------------------------\n");
    printf("Number of structures per space group:         %d \n", set.num_structures);
    printf("Number of molecules in the cell:              %d\n", set.Z);
    printf("Mean volume of unit cell:                     %f\n", set.vol_mean);
    printf("Standard deviation of unit cell volume:       %f\n", set.vol_std);
    printf("Spacegroup distribution type:                 %s\n", set.spg_dist_type);
    //printf("Specific radius proportion:                   %f\n", *sr );
    printf("Maximum attempts per space group:             %ld\n",set.max_attempts);
    printf("Volume attempts:                              %d\n", set.vol_attempts);
    printf("Random seed:                                  %d\n", set.random_seed);
    printf("Lattice angle standard deviation:             %f\n", set.angle_std);
    printf("Lattice principal component deviation:        %f\n", set.norm_dev);
    printf("Tolerance:                                    %f\n", set.tol);
    printf("Use rigid press:                              %d\n", set.rigid_press);
    if(set.n_mol_types > 1)
    {
        printf("Number of Molecules:                          %d\n", set.n_mol_types);
        printf("Stochiometry:                                 ");
        for(int i = 0; i < set.n_mol_types; i++)
        {
            printf("%d", set.stoic[i]);
            if(i != set.n_mol_types - 1)
                printf(":");
        }
        printf("\n");
    }
    printf("-----------------------------\n\n");

}

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
                          	int *random_seed)
{
    *Zp_max = 192; //useless argument
    printf("INPUT SETTINGS:\n");
    printf("-----------------------------\n");
    printf("Number of structures per layer group:         %d \n", *num_structures);
    printf("Number of molecules in the cell:              %d\n", *Z);
    printf("Mean volume of unit cell:                     %f\n", *volume_mean);
    printf("Standard deviation of unit cell volume:       %f\n", *volume_std);
    printf("Mean interface area:                          %f\n", *interface_area_mean);
    printf("Standard deviation of interface area:         %f\n", *interface_area_std);
    printf("Volume Multiplier:                            %d\n", *volume_multiplier);
    printf("Layergroup distribution type:                 %s\n", spg_dist_type);
    printf("Lattice vector:                               [[%f,%f,%f],[%f,%f,%f]]\n",global_lattice_vector_2d[0][0],
                                                global_lattice_vector_2d[0][1],global_lattice_vector_2d[0][2],
                                                global_lattice_vector_2d[1][0],global_lattice_vector_2d[1][1],
                                                global_lattice_vector_2d[1][2]  );
    printf("Maximum attempts per layer group:             %ld\n", *max_attempts);
    printf("Volume attempts:                              %d\n", *vol_attempt);
    printf("Random seed:                                  %d\n", *random_seed);
    printf("Tolerance:                                    %f\n", TOL);
    printf("-----------------------------\n\n");

}


void print_input_geometry(molecule* mol)
{
    printf("MOLECULAR GEOMETRY:\n");
    printf("-----------------------------\n");
    print_molecule(mol);
    printf("-----------------------------\n");
    printf("\n");
}


