#ifndef _ASU_UTILS_H_
#define _ASU_UTILS_H_
#include <stdio.h>
#include "asu.h"

// asu = asymmetric unit
void asu_init(asu *unit, int *stoic, int *n_atoms_in_mol, int n_mol_types);
void asu_allocate(asu *unit, int total_atoms);
void print_asu(asu *unit);
void asu_free(asu *unit);
void write_asu_extxyz(asu *unit, FILE *out, int structure_number);
void recenter_asu(asu *unit);

#endif
