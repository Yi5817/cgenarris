#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "algebra.h"
#include "check_structure.h"
#include "crystal.h"
#include "crystal_utils.h"
#include "molecule.h"
#include "read_input.h"
#include "pygenarris_mpi.h"
#include "spglib.h"

enum status{FAIL, SUCCESS};

#define MIN_STRUCTURES 100     
#define VOLUME_SIGMAS  5.0     // |det(lattice)| within mean +/- this many std
#define RIGID_TOL      0.02f   // Angstrom, intramolecular distances
#define SYMM_TOL       1e-3f   // fractional, symmetry images of atoms

// Declarations
int regression_test(MPI_Comm world_comm);
int check_geometry_out(molecule *mol, float *vdw_matrix, int dim, int Z,
                       double volume_mean, double volume_std);
static int check_composition(crystal *xtal, molecule *mol, int Z);
static int check_volume(crystal *xtal, double volume_mean, double volume_std);
static int check_rigid_molecules(crystal *xtal, molecule *mol);
static int check_symmetry(crystal *xtal);
char* get_metadata_line(char* line);
int read_geometry_out(FILE *fptr, crystal *xtal);

int regression_test(MPI_Comm world_comm)
{
    printf("Running regression test...\n");
    // settings
    molecule mol;
    float sr = 0.85;
    double volume_mean = 600;
    double volume_std = 40;
    double tol = 0.1;
    int num_structures = 2;
    int Z = 2;
    int seed = 223344;
    int volume_attempts = 10000;
    long max_attempts = 100000;
    char spg_dist_type[10] = "uniform";
    float norm_dev = 0.4;
    float angle_std = 8;

    // Read molecule
    read_geometry(&mol, "geometry.in");

    // Create vdw distance cutoff matrix
    int num_atoms_in_molecule = mol.num_of_atoms;
    int dim_vdw_matrix = num_atoms_in_molecule * Z ;
    float *vdw_cutoff_matrix = (float *) malloc( dim_vdw_matrix *
                                dim_vdw_matrix *
                                sizeof(float) ); //square matrix
    create_vdw_matrix_from_sr(&mol, vdw_cutoff_matrix, sr, Z);

    //call the generator from pygenarris_mpi
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
        volume_attempts,
        seed,
        norm_dev,
        angle_std,
        world_comm);

    int status = check_geometry_out(&mol, vdw_cutoff_matrix, dim_vdw_matrix,
                                    Z, volume_mean, volume_std);
    free(vdw_cutoff_matrix);
    return status;
}

int check_geometry_out(molecule *mol, float *vdw_matrix, int dim, int Z,
                       double volume_mean, double volume_std)
{
    crystal xtal;
    FILE *fptr = fopen("geometry.out", "r");
    if(!fptr)
    {
        printf("***ERROR: Cannot open geometry.out\n");
        return FAIL;
    }

    int n = 0;
    int n_bad = 0;
    while(read_geometry_out(fptr, &xtal))
    {
        n++;
        int ok = check_composition(&xtal, mol, Z);
        if(ok)
        {
            ok &= check_volume(&xtal, volume_mean, volume_std);
            ok &= check_rigid_molecules(&xtal, mol);
            ok &= check_symmetry(&xtal);
            if(!check_structure_with_vdw_matrix(xtal, vdw_matrix, dim, dim))
            {
                printf("FAIL structure %d (spg %d): vdw check rejects it\n",
                       n, xtal.spg);
                ok = 0;
            }
        }
        if(!ok)
            n_bad++;

        free(xtal.Xcord);
        free(xtal.Ycord);
        free(xtal.Zcord);
        free(xtal.atoms);
    }
    fclose(fptr);

    printf("Checked %d structures, %d failed\n", n, n_bad);
    if(n < MIN_STRUCTURES)
    {
        printf("FAIL: expected at least %d structures\n", MIN_STRUCTURES);
        return FAIL;
    }
    return n_bad ? FAIL : SUCCESS;
}

// Z and the molecule size match the request, and the atoms are the input
// molecule repeated Z times in order.
static int check_composition(crystal *xtal, molecule *mol, int Z)
{
    int N = mol->num_of_atoms;
    if(xtal->Z != Z || xtal->num_atoms_in_molecule != N)
    {
        printf("FAIL structure spg %d: Z=%d napm=%d, expected Z=%d napm=%d\n",
               xtal->spg, xtal->Z, xtal->num_atoms_in_molecule, Z, N);
        return 0;
    }
    for(int i = 0; i < Z*N; i++)
    {
        if(xtal->atoms[2*i] != mol->atoms[2*(i%N)] ||
           xtal->atoms[2*i+1] != mol->atoms[2*(i%N)+1])
        {
            printf("FAIL structure spg %d: atom %d is %c%c, expected %c%c\n",
                   xtal->spg, i, xtal->atoms[2*i], xtal->atoms[2*i+1],
                   mol->atoms[2*(i%N)], mol->atoms[2*(i%N)+1]);
            return 0;
        }
    }
    return 1;
}

// Cell volume is finite and within the sampled volume distribution.
static int check_volume(crystal *xtal, double volume_mean, double volume_std)
{
    double vol = fabs(det_mat3b3(xtal->lattice_vectors));
    double lo = volume_mean - VOLUME_SIGMAS * volume_std;
    double hi = volume_mean + VOLUME_SIGMAS * volume_std;
    if(!isfinite(vol) || vol < lo || vol > hi)
    {
        printf("FAIL structure spg %d: volume %f outside [%f, %f]\n",
               xtal->spg, vol, lo, hi);
        return 0;
    }
    return 1;
}

// Every molecule in the cell is a rigid copy of the input molecule:
// all intramolecular distances match geometry.in.
static int check_rigid_molecules(crystal *xtal, molecule *mol)
{
    int N = mol->num_of_atoms;
    for(int m = 0; m < xtal->Z; m++)
    {
        for(int i = 0; i < N; i++)
        for(int j = i+1; j < N; j++)
        {
            float ref[3] = {mol->X[i] - mol->X[j],
                            mol->Y[i] - mol->Y[j],
                            mol->Z[i] - mol->Z[j]};
            int a = m*N + i;
            int b = m*N + j;
            float got[3] = {xtal->Xcord[a] - xtal->Xcord[b],
                            xtal->Ycord[a] - xtal->Ycord[b],
                            xtal->Zcord[a] - xtal->Zcord[b]};
            float d = vector3_norm(got) - vector3_norm(ref);
            if(!isfinite(d) || fabsf(d) > RIGID_TOL)
            {
                printf("FAIL structure spg %d: molecule %d distorted, "
                       "atoms %d-%d off by %f A\n", xtal->spg, m, i, j, d);
                return 0;
            }
        }
    }
    return 1;
}

// The structure has the symmetry of the attempted space group: every
// operation maps each atom onto an atom of the same species (mod lattice).
static int check_symmetry(crystal *xtal)
{
    int rotations[192][3][3];
    double translations[192][3];
    int hall_number = hall_number_from_spg(xtal->spg);
    int num_ops = spg_get_symmetry_from_database(rotations, translations,
                                                 hall_number);
    int total_atoms = xtal->Z * xtal->num_atoms_in_molecule;
    int ok = 1;

    convert_xtal_to_fractional(xtal);
    for(int op = 0; op < num_ops && ok; op++)
    {
        float trans[3] = {translations[op][0], translations[op][1],
                          translations[op][2]};
        for(int a = 0; a < total_atoms && ok; a++)
        {
            float image[3] = {xtal->Xcord[a], xtal->Ycord[a], xtal->Zcord[a]};
            vector3_intmat3b3_multiply(rotations[op], image, image);
            vector3_add(trans, image, image);

            int found = 0;
            for(int b = 0; b < total_atoms && !found; b++)
            {
                if(xtal->atoms[2*b] != xtal->atoms[2*a] ||
                   xtal->atoms[2*b+1] != xtal->atoms[2*a+1])
                    continue;
                float diff[3] = {image[0] - xtal->Xcord[b],
                                 image[1] - xtal->Ycord[b],
                                 image[2] - xtal->Zcord[b]};
                found = fabsf(diff[0] - roundf(diff[0])) < SYMM_TOL &&
                        fabsf(diff[1] - roundf(diff[1])) < SYMM_TOL &&
                        fabsf(diff[2] - roundf(diff[2])) < SYMM_TOL;
            }
            if(!found)
            {
                printf("FAIL structure spg %d: symmetry op %d/%d has no "
                       "image for atom %d\n", xtal->spg, op, num_ops, a);
                ok = 0;
            }
        }
    }
    convert_xtal_to_cartesian(xtal);
    return ok;
}

int main(int argc, char **argv)
{
    //Initialise MPI
    MPI_Init(&argc, &argv);
    MPI_Comm world_comm = MPI_COMM_WORLD;

    int status;
    status = regression_test(world_comm);
    if(status == SUCCESS)
    {
        printf("Regression tests passed\n");
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    else if(status == FAIL)
    {
        printf("Regression test failed\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    MPI_Finalize();
    exit(EXIT_FAILURE);
}


int read_geometry_out(FILE* fptr, crystal *xtal)
{
  char *line = NULL;
  char *sub_line = NULL;
  size_t len = 0;
  int num_atom = 0;

  int spg = 0;
  int Z = 0;
  int wyckoff_position = 0;
  int Zp = 0;
  int i = 0;
  int counter = 0;
  int num_atoms_in_molecule = 0;

  if(!fptr)
  {
       printf("***ERROR: Cannot open the file \n");
       exit(EXIT_FAILURE);
  }

    while(getline(&line,&len,fptr) != -1)
    {
        // Starts with # then
        if(strstr(line,"#") != NULL || *line== '\n')
        {
            if(strstr(line, "attempted_spacegroup") != NULL)
            {
                 sub_line = get_metadata_line(line);
                 spg = atoi(sub_line);
                 // Set temporary counting index for new xtal structure
                 i = 0;
                 counter = 0;
            }

            else if (strstr(line, "Z") != NULL)
            {
                sub_line = get_metadata_line(line);
                Z = atoi(sub_line);
            }

            else if (strstr(line,"number_of_atoms_in_molecule") != NULL)
            {
                sub_line = get_metadata_line(line);
                num_atoms_in_molecule = atof(sub_line);
            }

            else if (strstr(line,"END") != NULL)
            {
                 xtal->spg = spg;
                 xtal->wyckoff_position = wyckoff_position;
                 xtal->num_atoms_in_molecule = num_atoms_in_molecule;
                 xtal->Z = Z;
                 xtal->Zp = Zp;
                 return 1;
            }

            continue;
        }

        sub_line = strtok(line, "    ");
        if(strcmp(sub_line, "lattice_vector") == 0)
        {
            sub_line = strtok(NULL, "    ");
            xtal->lattice_vectors[i][0] = atof(sub_line);
            sub_line = strtok(NULL, "    ");
            xtal->lattice_vectors[i][1] = atof(sub_line);
            sub_line = strtok(NULL, "    ");
            xtal->lattice_vectors[i][2] = atof(sub_line);
            if(!i)
            {
                num_atom = num_atoms_in_molecule *Z;
                xtal->Xcord = (float *) malloc (num_atom*sizeof(float));
                xtal->Ycord = (float *) malloc (num_atom*sizeof(float));
                xtal->Zcord = (float *) malloc (num_atom*sizeof(float));
                xtal->atoms = (char * ) malloc (2*num_atom*sizeof(char));
             }
            i++;
        }

        if(strcmp(sub_line, "atom") == 0)
        {
            xtal->Xcord[counter] = atof(strtok(NULL, "    "));
            xtal->Ycord[counter] = atof(strtok(NULL, "    "));
            xtal->Zcord[counter] = atof(strtok(NULL, "    "));

            sub_line = strtok(NULL, " ");
            xtal->atoms[2*counter] = *sub_line;
            if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
                xtal->atoms[2*counter+1] = ' ';
            else
                xtal->atoms[2*counter+1] = *(sub_line+1);
            counter++;
        }
    }

    return 0;
}

// extract the meta data line like Z, spg etc
char* get_metadata_line(char* line)
{
    char* sub_line = NULL;
    sub_line = strtok(line," ");
    sub_line = strtok(NULL," ");
    sub_line = strtok(NULL," ");
    return sub_line;
}