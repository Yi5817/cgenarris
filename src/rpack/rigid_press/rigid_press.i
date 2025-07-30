%module rigid_press
%{
#include "rigid_press.h"
%}

%include "../../crystal.h"

%{
#define SWIG_FILE_WITH_INIT
%}
%include "../../numpy.i"
%init %{
import_array();
%}


//Settings for optmization
typedef struct
{
    int max_iteration;
    int cell_family;   //See keys above
    int spg;
    
} Opt_settings;

typedef enum {
    SUCCESS,           // All good
    ITER_LIMIT,        // Reached max iteration limit 
    FALSE_CONVERGENCE, // Optimization ended, but gradient is still large
    MISC_FAILURE       // Some other reason for failure
} Opt_status;

%apply (float* IN_ARRAY1, int DIM1) {(float *cutoff_matrix, int placeholder)};
// optimizes a molecular crystal using a regularized rigid-body interaction
Opt_status optimize_crystal(crystal *xtl, // a molecular crystal in the Genarris crystal format [1]
			    float *cutoff_matrix, // distance cutoff between pairs of atoms in the crystallized molecule [(xtl->Z*xtl->num_atoms_in_molecule)^2]
			    int placeholder, // Used for swig interface
			    Opt_settings set);
