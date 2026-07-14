#include "asu_utils.h"
#include "asu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "algebra.h"



void asu_init(asu *unit, int *stoic, int *n_atoms_in_mol, int n_mol_types)
{

    int n_atoms = 0;
    int n_mol_asym = 0;  // Number of molecules in asymmetric unit

    // Find number of atoms first
    for(int m = 0; m < n_mol_types; m++)
    {
        // Number of atoms in asym unit
        n_atoms += n_atoms_in_mol[m] * stoic[m];
        n_mol_asym += stoic[m];
    }


    // Allocate memory
    int tbytes = n_mol_types * sizeof(int);  // Total number of mol types
    int mbytes = n_mol_asym  * sizeof(int);  // Total num of molecules
    unit->stoic            = (int *) malloc(tbytes);
    unit->n_atoms_in_mol   = (int *) malloc(tbytes);
    unit->mol_index = (int *) malloc(mbytes);
    unit->mol_types = (int *) malloc(mbytes);
    unit->n_mols = n_mol_asym;
    asu_allocate(unit, n_atoms);

    // store details
    memcpy(unit->stoic, stoic, tbytes);
    memcpy(unit->n_atoms_in_mol, n_atoms_in_mol, tbytes);
    unit->n_mol_types = n_mol_types;
    unit->n_atoms = n_atoms;

    // Get molecule index
    int at = 0;
    int mol_id = 0;

    for(int m = 0; m < n_mol_types; m++)
       {
       for(int st = 0; st < unit->stoic[m]; st++)
              {
                     unit->mol_types[mol_id] = m;
                     unit->mol_index[mol_id] = at;
                     at += n_atoms_in_mol[m];
                     mol_id++;
              }
       }

}


void asu_allocate(asu *unit, int total_atoms)
{
    int fbytes = total_atoms * sizeof(float);
    unit->Xcord = malloc(fbytes);
    unit->Ycord = malloc(fbytes);
    unit->Zcord = malloc(fbytes);
    unit->com   = malloc(3 * unit->n_mols * sizeof(float));
    unit->atoms = malloc(total_atoms * sizeof(char) *2);
}

void print_asu(asu *unit)
{
	int N = (*unit).n_atoms;
	printf("#total atoms in the asymmetric unit = %d \n\n", N);

	for(int i =0; i < N; i++)
	{
		printf("atom %12f %12f %12f %4c%c \n",
			(*unit).Xcord[i], (*unit).Ycord[i], (*unit).Zcord[i], (*unit).atoms[2*i],(*unit).atoms[2*i+1]);
	}
}

void asu_free(asu *unit)
{
    free(unit->Xcord);
    free(unit->Ycord);
    free(unit->Zcord);
    free(unit->atoms);
    free(unit->com);
    free(unit->stoic);
    free(unit->n_atoms_in_mol);
    free(unit->mol_types);
    free(unit->mol_index);

}


/*
Writes one asymmetric unit as a single extended-XYZ block.

Args:
    unit: The asymmetric unit to serialize.
    out: Open, writable file handle (blocks are appended in order).
    structure_number: 1-based index of this structure within the file.

The block has the standard extxyz layout: an atom-count line, a
comment line of space-separated key=value metadata, and one
"element x y z" row per atom. No Lattice key is written and pbc is
"F F F", so ase.io.read loads each block as a non-periodic Atoms
object with the metadata available in Atoms.info.
*/
void write_asu_extxyz(asu *unit, FILE *out, int structure_number)
{
    fprintf(out, "%d\n", unit->n_atoms);

    fprintf(out, "Properties=species:S:1:pos:R:3 pbc=\"F F F\"");
    fprintf(out, " structure_number=%d", structure_number);
    fprintf(out, " n_atoms=%d n_mols=%d n_mol_types=%d",
            unit->n_atoms, unit->n_mols, unit->n_mol_types);
    fprintf(out, " sr=%.6f", unit->sr);

    fprintf(out, " stoic=\"");
    for(int i = 0; i < unit->n_mol_types; i++)
        fprintf(out, "%d%s", unit->stoic[i], i + 1 < unit->n_mol_types ? " " : "");
    fprintf(out, "\"");

    fprintf(out, " mol_types=\"");
    for(int i = 0; i < unit->n_mols; i++)
        fprintf(out, "%d%s", unit->mol_types[i], i + 1 < unit->n_mols ? " " : "");
    fprintf(out, "\"\n");

    for(int at = 0; at < unit->n_atoms; at++)
    {
        char e0 = unit->atoms[2*at];
        char e1 = unit->atoms[2*at + 1];
        if(e1 == ' ' || e1 == '\0' || e1 == '\n')
            fprintf(out, "%c %.8f %.8f %.8f\n",
                    e0, unit->Xcord[at], unit->Ycord[at], unit->Zcord[at]);
        else
            fprintf(out, "%c%c %.8f %.8f %.8f\n",
                    e0, e1, unit->Xcord[at], unit->Ycord[at], unit->Zcord[at]);
    }
    fflush(out);
}

void recenter_asu(asu *unit){

    int N = (*unit).n_atoms;
	float xcom = 0, ycom = 0 , zcom = 0;
	for(int i = 0; i < N; i++)
	{
		xcom += (*unit).Xcord[i];
		ycom += (*unit).Ycord[i];
		zcom += (*unit).Zcord[i];
	}

	xcom /= N;
	ycom /= N;
	zcom /= N;

	for(int i = 0; i < N; i++)
	{
		(*unit).Xcord[i] -= xcom;
		(*unit).Ycord[i] -= ycom;
		(*unit).Zcord[i] -= zcom;
	}

}
