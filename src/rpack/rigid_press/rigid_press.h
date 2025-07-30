// Interface for the fast, approximate molecular crystal optimizer

#ifndef RIGID_PRESS_H
#define RIGID_PRESS_H

#include "../../crystal.h"

// crystal families & their lattice vector constraints
#define TRICLINIC    0 // no constraints on lattice vectors
#define MONOCLINIC   1 // (a, 0, 0), (0, b, 0), (c, 0, d)
#define ORTHORHOMBIC 2 // (a, 0, 0), (0, b, 0), (0, 0, c)
#define TETRAGONAL   3 // (a, 0, 0), (0, a, 0), (0, 0, b)
#define HEXAGONAL    4 // (a, 0, 0), (-a/2, a*sqrt(3)/2, 0), (0, 0, b)
#define CUBIC        5 // (a, 0, 0), (0, a, 0), (0, 0, a)

//Settings for optmization
typedef struct
{
    int max_iteration;
    int cell_family;   //See keys above
    int spg;
    
} Opt_settings;

typedef enum Opt_status
{
    SUCCESS,           // All good
    ITER_LIMIT,        // Reached max iteration limit 
    FALSE_CONVERGENCE, // Optimization ended, but gradient is still large
    MISC_FAILURE       // Some other reason for failure
} Opt_status;

// optimizes a molecular crystal using a regularized rigid-body interaction
Opt_status optimize_crystal(crystal *xtl, // a molecular crystal in the Genarris crystal format [1]
			    float *cutoff_matrix, // distance cutoff between pairs of atoms in the crystallized molecule [(xtl->Z*xtl->num_atoms_in_molecule)^2]
			    int placeholder,
			    Opt_settings set);

#endif
