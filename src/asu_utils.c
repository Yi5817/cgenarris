#include "asu_utils.h"
#include "cgenarris/asu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "check_structure.h"

int asu_init(asu *unit, const molecule *mol, const int *stoic, int n_mol_types)
{
    memset(unit, 0, sizeof(*unit));

    if(n_mol_types < 1)
    {
        fprintf(stderr, "***ERROR: asu_init: n_mol_types must be >= 1, got %d\n",
                n_mol_types);
        return -1;
    }

    int n_atoms = 0;
    int n_mols = 0;
    for(int m = 0; m < n_mol_types; m++)
    {
        if(stoic[m] < 1)
        {
            fprintf(stderr, "***ERROR: asu_init: stoichiometry[%d] must be >= 1, got %d\n",
                    m, stoic[m]);
            return -1;
        }
        if(mol[m].num_of_atoms < 1)
        {
            fprintf(stderr, "***ERROR: asu_init: molecule %d has %d atoms\n",
                    m, mol[m].num_of_atoms);
            return -1;
        }
        n_atoms += mol[m].num_of_atoms * stoic[m];
        n_mols += stoic[m];
    }

    if(n_mols < 2)
    {
        fprintf(stderr, "***ERROR: asu_init: an asymmetric unit needs at least "
                "two molecules (sum of stoichiometry), got %d\n", n_mols);
        return -1;
    }

    unit->n_atoms = n_atoms;
    unit->n_mols = n_mols;
    unit->n_mol_types = n_mol_types;
    unit->Xcord = (float *)malloc(n_atoms * sizeof(float));
    unit->Ycord = (float *)malloc(n_atoms * sizeof(float));
    unit->Zcord = (float *)malloc(n_atoms * sizeof(float));
    unit->vdw_radii = (float *)malloc(n_atoms * sizeof(float));
    unit->atoms = (char *)malloc(2 * n_atoms * sizeof(char));
    unit->mol_index = (int *)malloc(n_mols * sizeof(int));
    unit->mol_types = (int *)malloc(n_mols * sizeof(int));
    unit->n_atoms_in_mol = (int *)malloc(n_mol_types * sizeof(int));
    unit->stoic = (int *)malloc(n_mol_types * sizeof(int));

    // Molecule tables and species, ordered by type then copy.
    int at = 0;
    int k = 0;
    for(int m = 0; m < n_mol_types; m++)
    {
        unit->stoic[m] = stoic[m];
        unit->n_atoms_in_mol[m] = mol[m].num_of_atoms;
        for(int copy = 0; copy < stoic[m]; copy++)
        {
            unit->mol_types[k] = m;
            unit->mol_index[k] = at;
            memcpy(unit->atoms + 2 * at, mol[m].atoms,
                   2 * mol[m].num_of_atoms * sizeof(char));
            at += mol[m].num_of_atoms;
            k++;
        }
    }

    for(int i = 0; i < n_atoms; i++)
    {
        if(atom_vdw_radius(unit->atoms[2 * i], unit->atoms[2 * i + 1],
                           unit->vdw_radii + i))
        {
            fprintf(stderr, "***ERROR: asu_init: no van der Waals radius for "
                    "species '%c%c' (atom %d)\n",
                    unit->atoms[2 * i], unit->atoms[2 * i + 1], i);
            asu_free(unit);
            return -1;
        }
    }

    return 0;
}

void asu_free(asu *unit)
{
    free(unit->Xcord);
    free(unit->Ycord);
    free(unit->Zcord);
    free(unit->vdw_radii);
    free(unit->atoms);
    free(unit->mol_index);
    free(unit->mol_types);
    free(unit->n_atoms_in_mol);
    free(unit->stoic);
    memset(unit, 0, sizeof(*unit));
}

void asu_recenter(asu *unit)
{
    int N = unit->n_atoms;
    float xcom = 0, ycom = 0, zcom = 0;
    for(int i = 0; i < N; i++)
    {
        xcom += unit->Xcord[i];
        ycom += unit->Ycord[i];
        zcom += unit->Zcord[i];
    }
    xcom /= N;
    ycom /= N;
    zcom /= N;
    for(int i = 0; i < N; i++)
    {
        unit->Xcord[i] -= xcom;
        unit->Ycord[i] -= ycom;
        unit->Zcord[i] -= zcom;
    }
}

void asu_print(const asu *unit)
{
    printf("#total atoms in the asymmetric unit = %d\n\n", unit->n_atoms);
    for(int i = 0; i < unit->n_atoms; i++)
    {
        printf("atom %12f %12f %12f %c%c\n",
               unit->Xcord[i], unit->Ycord[i], unit->Zcord[i],
               unit->atoms[2 * i], unit->atoms[2 * i + 1]);
    }
}

// Prints "#key = v0 v1 ..." for an int array.
static void write_int_list(FILE *out, const char *key, const int *values, int n)
{
    fprintf(out, "#%s =", key);
    for(int i = 0; i < n; i++)
        fprintf(out, " %d", values[i]);
    fprintf(out, "\n");
}

void asu_write_block(const asu *unit, FILE *out, int structure_number)
{
    fprintf(out, "####### BEGIN STRUCTURE #######\n");
    fprintf(out, "#structure_number = %d\n", structure_number);
    fprintf(out, "#number_of_atoms = %d\n", unit->n_atoms);
    fprintf(out, "#number_of_molecules = %d\n", unit->n_mols);
    fprintf(out, "#number_of_molecule_types = %d\n", unit->n_mol_types);
    write_int_list(out, "stoichiometry", unit->stoic, unit->n_mol_types);
    write_int_list(out, "number_of_atoms_in_molecule_type", unit->n_atoms_in_mol,
                   unit->n_mol_types);
    write_int_list(out, "molecule_types", unit->mol_types, unit->n_mols);
    write_int_list(out, "molecule_index", unit->mol_index, unit->n_mols);
    fprintf(out, "#sr = %f\n", unit->sr);
    fprintf(out, "#\"All distances in Angstroms and using Cartesian coordinate system\"\n");

    for(int i = 0; i < unit->n_atoms; i++)
    {
        fprintf(out, "atom %12f %12f %12f  %c%c \n", unit->Xcord[i],
                unit->Ycord[i], unit->Zcord[i],
                unit->atoms[2 * i], unit->atoms[2 * i + 1]);
    }
    fprintf(out, "#######  END  STRUCTURE #######\n\n");
    fflush(out);
}
