/*
 * Standalone test for check_structure_with_vdw_matrix.
 *
 * Prints PASS / FAIL and the verdict integer.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "crystal.h"
#include "molecule.h"
#include "read_input.h"
#include "check_structure.h"

static int parse_structure_file(const char *path, crystal *xtal);
static float *read_vdw_matrix(const char *path, int *dim_out);

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        fprintf(stderr,
            "Usage: %s <geometry.in> <vdw_matrix.txt> <geometry.out>\n",
            argv[0]);
        return 2;
    }

    const char *mol_path = argv[1];
    const char *vdw_path = argv[2];
    const char *xtal_path = argv[3];

    molecule mol;
    read_geometry(&mol, (char *)mol_path);
    printf("Molecule: %d atoms read from %s\n", mol.num_of_atoms, mol_path);

    int dim = 0;
    float *vdw_matrix = read_vdw_matrix(vdw_path, &dim);
    printf("vdw_matrix: %d x %d read from %s\n", dim, dim, vdw_path);

    crystal xtal;
    if (!parse_structure_file(xtal_path, &xtal))
    {
        fprintf(stderr, "***ERROR: failed to parse %s\n", xtal_path);
        return 2;
    }
    xtal.num_atoms_in_molecule = mol.num_of_atoms;
    int total_atoms = xtal.Z * xtal.num_atoms_in_molecule;
    printf("Crystal: Z=%d, atoms_per_mol=%d, total_atoms=%d\n",
        xtal.Z, xtal.num_atoms_in_molecule, total_atoms);
    printf("vdw_matrix dim vs total_atoms: %d vs %d %s\n",
        dim, total_atoms, dim == total_atoms ? "(match)" : "(MISMATCH)");

    int verdict = check_structure_with_vdw_matrix(xtal, vdw_matrix, dim, dim);
    printf("\nverdict = %d  ->  %s\n", verdict,
        verdict ? "PASS (structure accepted)" : "FAIL (structure rejected)");

    free(vdw_matrix);
    free(xtal.Xcord);
    free(xtal.Ycord);
    free(xtal.Zcord);
    free(xtal.atoms);
    free(mol.X);
    free(mol.Y);
    free(mol.Z);
    free(mol.atoms);

    return verdict ? 0 : 1;
}

static float *read_vdw_matrix(const char *path, int *dim_out)
{
    FILE *fp = fopen(path, "r");
    if (!fp)
    {
        fprintf(stderr, "***ERROR: cannot open %s\n", path);
        exit(2);
    }

    // Two-pass: count floats, then read into a square matrix
    int count = 0;
    float dummy;
    while (fscanf(fp, "%f", &dummy) == 1)
        count++;

    int dim = (int)(sqrtf((float)count) + 0.5f);
    if (dim * dim != count || dim <= 0)
    {
        fprintf(stderr,
            "***ERROR: %s has %d floats; not a perfect square\n",
            path, count);
        exit(2);
    }

    rewind(fp);
    float *m = (float *)malloc(sizeof(float) * dim * dim);
    for (int i = 0; i < dim * dim; i++)
    {
        if (fscanf(fp, "%f", &m[i]) != 1)
        {
            fprintf(stderr,
                "***ERROR: %s second pass failed at index %d\n", path, i);
            exit(2);
        }
    }
    fclose(fp);
    *dim_out = dim;
    return m;
}

/*
 * Minimal parser for a single structure block in the cgenarris geometry.out
 * format. Reads Z, num_atoms_in_molecule, lattice_vectors, and atom lines.
 * Stops at "END STRUCTURE" or EOF.
 */
static int parse_structure_file(const char *path, crystal *xtal)
{
    FILE *fp = fopen(path, "r");
    if (!fp)
    {
        fprintf(stderr, "***ERROR: cannot open %s\n", path);
        return 0;
    }

    char *line = NULL;
    size_t len = 0;
    int lat_row = 0;
    int atom_count = 0;
    int allocated = 0;
    int Z = 0;
    int napm = 0;
    int found_begin = 0;

    while (getline(&line, &len, fp) != -1)
    {
        if (strstr(line, "BEGIN STRUCTURE") != NULL)
        {
            found_begin = 1;
            continue;
        }
        if (strstr(line, "END STRUCTURE") != NULL)
            break;

        // Metadata lines start with '#'
        if (line[0] == '#')
        {
            char *eq = strchr(line, '=');
            if (!eq) continue;
            eq++;
            while (*eq == ' ') eq++;
            if (strstr(line, "Z =") != NULL ||
                strstr(line, "Z=") != NULL)
            {
                // First "Z =" line is Z. (Other 'Z'-containing keys like
                // "SPGLIB" are handled by the strstr above being narrow.)
                if (Z == 0) Z = atoi(eq);
            }
            else if (strstr(line, "number_of_atoms_in_molecule") != NULL)
            {
                napm = atoi(eq);
            }
            continue;
        }

        // Allocate once we know dimensions
        if (!allocated && Z > 0 && napm > 0)
        {
            int total = Z * napm;
            xtal->Xcord = (float *)malloc(sizeof(float) * total);
            xtal->Ycord = (float *)malloc(sizeof(float) * total);
            xtal->Zcord = (float *)malloc(sizeof(float) * total);
            xtal->atoms = (char *)malloc(sizeof(char) * 2 * total);
            allocated = 1;
        }

        char *tok = strtok(line, " \t\n");
        if (!tok) continue;

        if (strcmp(tok, "lattice_vector") == 0 && lat_row < 3)
        {
            xtal->lattice_vectors[lat_row][0] = atof(strtok(NULL, " \t\n"));
            xtal->lattice_vectors[lat_row][1] = atof(strtok(NULL, " \t\n"));
            xtal->lattice_vectors[lat_row][2] = atof(strtok(NULL, " \t\n"));
            lat_row++;
        }
        else if (strcmp(tok, "atom") == 0 && allocated)
        {
            xtal->Xcord[atom_count] = atof(strtok(NULL, " \t\n"));
            xtal->Ycord[atom_count] = atof(strtok(NULL, " \t\n"));
            xtal->Zcord[atom_count] = atof(strtok(NULL, " \t\n"));
            char *sym = strtok(NULL, " \t\n");
            xtal->atoms[2 * atom_count] = sym ? sym[0] : ' ';
            xtal->atoms[2 * atom_count + 1] =
                (sym && sym[1] && sym[1] != '\n') ? sym[1] : ' ';
            atom_count++;
        }
    }

    free(line);
    fclose(fp);

    if (!found_begin || Z == 0 || napm == 0)
        return 0;

    xtal->Z = Z;
    xtal->num_atoms_in_molecule = napm;
    return 1;
}
