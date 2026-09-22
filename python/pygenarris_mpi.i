%module pygenarris_mpi
%{
#include "cgenarris/read_input.h"
#include "spg_generation.h"
#include "pygenarris.h"
#include "cgenarris/pygenarris_mpi.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "randomgen.h"
#include "spglib.h"
#include "mpi.h"
#include "cgenarris_mpi.h"
#include "cgenarris/asu_generation.h"

%}

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}

%include "mpi4py/mpi4py.i"
%mpi4py_typemap(Comm, MPI_Comm);

%apply ( float* IN_ARRAY2, int DIM1, int DIM2) {(float *vdw_matrix, int dim1, int dim2)};

void mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
    float *vdw_matrix,
    int dim1,
    int dim2,
    int num_structures,
    int Z,
    double volume_mean1,
    double volume_std1,
    double tol1,
    long max_attempts,
    char *spg_dist_type,
    int vol_attempt,
    int random_seed,
    float norm_dev,
    float angle_std,
    MPI_Comm world_comm);

%apply (float INPLACE_ARRAY2[ANY][ANY]) {(float lattice_vector_2d_from_geo[2][3])};
void mpi_generate_layer_with_vdw_cutoff_matrix(
	float *vdw_matrix,
	int dim1,
	int dim2,
	int num_structures,
	int Z,
	double volume_mean1,
	double volume_std1,
	double interface_area_mean1,
	double interface_area_std1,
	int volume_multiplier,
	double tol1,
	long max_attempts,
	char *spg_dist_type,
	float lattice_vector_2d_from_geo[2][3],
        int vol_attempt,
        int random_seed,
	MPI_Comm world_comm);

void get_compatible_spg(
    int Z,
    const char *geometry_file);

int num_compatible_spacegroups(int Z, double tolerance);

%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *positions, int n_atoms_total, int ncols)};
%apply (int* IN_ARRAY1, int DIM1) {(int *n_atoms_per_mol, int n_mol_types)};
%apply (int* IN_ARRAY1, int DIM1) {(int *stoichiometry, int n_stoichiometry)};
%feature("autodoc", "generate_asymmetric_units(positions, species, n_atoms_per_mol, stoichiometry, "
    "num_structures, sr_min, sr_max, max_attempts, random_seed, output_file, comm) -> int\n\n"
    "Generate random asymmetric units in parallel and write them as geometry.out blocks.\n\n"
    "positions:       (n_atoms, 3) float64 array, all molecule types concatenated, Angstrom\n"
    "species:         str, two chars per atom, space padded ('C H Cl' -> 'C H Cl')\n"
    "n_atoms_per_mol: int32 array, atoms of each molecule type\n"
    "stoichiometry:   int32 array, copies of each molecule type per unit (all >= 1, sum >= 2)\n"
    "num_structures:  total number of asymmetric units\n"
    "sr_min, sr_max:  open window for the closest intermolecular d_ij/(r_i+r_j)\n"
    "max_attempts:    placements tried per asymmetric unit before giving up\n"
    "random_seed:     base seed; rank r uses seed + r; 0 = time based\n"
    "output_file:     path of the merged output file (overwritten)\n"
    "comm:            mpi4py communicator\n\n"
    "Returns the number of asymmetric units written, or -1 on invalid input.")
    generate_asymmetric_units;
int generate_asymmetric_units(
    double *positions,
    int n_atoms_total,
    int ncols,
    char *species,
    int *n_atoms_per_mol,
    int n_mol_types,
    int *stoichiometry,
    int n_stoichiometry,
    int num_structures,
    double sr_min,
    double sr_max,
    long max_attempts,
    int random_seed,
    char *output_file,
    MPI_Comm world_comm);

//void send_xtal(MPI_Comm comm, int destination, crystal* xtal, int total_atoms);
//void receive_xtal(MPI_Comm comm, int source, crystal* xtal, int total_atoms);
void find_allowed_positions_using_molecular_symmetry(char mol_sym[6],
	int Z, int Zpp);

void allocate_xtal(crystal* xtal, int Z, int N);

%include "cgenarris/crystal.h"
%apply (double INPLACE_ARRAY2[ANY][ANY]) {(double lattice_vector[3][3])};
%apply (double* IN_ARRAY1, int DIM1) {(double *Xc, int total_atoms1)};
%apply (double* IN_ARRAY1, int DIM1) {(double *Yc, int total_atoms2)};
%apply (double* IN_ARRAY1, int DIM1) {(double *Zc, int total_atoms3)};
void create_crystal_from_array(crystal *xtal, double lattice_vector[3][3],
			       double *Xc,int total_atoms1, double *Yc,int total_atoms2,
			       double *Zc, int total_atoms3, char *atoms,  int total_atoms,
			       int Z, int spg);
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *X, int total_atoms1)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *Y, int total_atoms2)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *Z, int total_atoms3)};
void create_array_from_crystal(crystal *xtal, double lattice_vector[3][3],
			       double *X,int total_atoms1, double *Y,int total_atoms2,
			       double *Z, int total_atoms3, char* atoms, int total_atoms,
			       int Z, int spg);

%apply (double INPLACE_ARRAY1[ANY]) {(double a[3])};

void print_crystal(crystal* xtal);
void free_xtal(crystal* xtal);


%apply ( float* IN_ARRAY2, int DIM1, int DIM2) {(float *vdw_matrix, int dim1, int dim2)};
int check_structure_with_vdw_matrix(crystal random_crystal,
	float *vdw_matrix,
	int dim1,
	int dim2);
