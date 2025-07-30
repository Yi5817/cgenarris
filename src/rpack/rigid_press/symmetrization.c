
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "../../spglib.h"
#include "d_algebra.h"
#include "../../crystal.h"
#include "symmetrization.h"

static void quaternion2matrix(double *quat, double *mat, const int inv);
static void matrix2quaternion(double *rot, double *quat, int *inv);
static void symmetrize_matrix(double *mat, int dim, int spg);
static void symmetrize_vector(double *vec, int dim, double lattice[3][3], int spg);

void compute_molecule_COM(crystal xtal, double com[3], int i);
void convert_xtal_to_cartesian(crystal *xtal);
void convert_xtal_to_fractional(crystal *xtal);

static void test_symmetrize_state();
static void symmetrize_gradient_pos(double *vec, int dim, double lattice[3][3], int spg);

static void symmetrize_gradient_pos(double *vec, int dim, double lattice[3][3], int spg)
{
    // Get symm operations
    double translations[192][3];
    int rotations[192][3][3];
    int hall_number = hall_number_from_spg(spg);
    int Z = spg_get_symmetry_from_database(rotations,
					   translations,
					   hall_number);
    assert(dim == Z);
    // Get inverse lattice
    double inverse_lattice[3][3], inverse_lattice_t[3][3], lattice_t[3][3];
    inverse_mat3b3(inverse_lattice, lattice);

    mat3b3_transpose(inverse_lattice_t, inverse_lattice);
    mat3b3_transpose(lattice_t, lattice);

    // Convert to fractional space
    double vec_frac[3*dim];
    for(int i = 0; i < dim; i++)
    {
	vector3_mat3b3_multiply(inverse_lattice_t, vec + 3*i, vec_frac +3*i);
    }

   
    // Array to compute symm average
    double vec_symm[3] = {0};
    
    for(int op = 0; op < Z; op++)
    {
	int (*rot_i)[3] = rotations[op];
	double rot[3][3] = {{rot_i[0][0], rot_i[0][1],rot_i[0][2]},
			   {rot_i[1][0], rot_i[1][1],rot_i[1][2]},
			   {rot_i[2][0], rot_i[2][1],rot_i[2][2]}};
	double trans[3] = {translations[op][0],
			  translations[op][1],
			  translations[op][2]};
	
	double inv_rot[3][3];
	inverse_mat3b3(inv_rot, rot);

	double temp[3];
	vector3_subtract(vec_frac + 3*op, trans, vec_frac + 3*op);
	vector3_mat3b3_multiply(inv_rot, vec_frac + 3*op, temp);
	vector3_add(temp, vec_symm, vec_symm);
    }

    // Take mean over total operations
    for(int i = 0; i < 3; i++)
	vec_symm[i] /= Z;
    
    // Apply symmetry operation to regenerate the entire vec
    for(int i = 0; i < Z; i++)
    {
	int (*rot)[3] = rotations[i];
	vector3_intmat3b3_multiply(rot, vec_symm, vec_frac + 3*i);
	double trans[3] = {translations[i][0],
			  translations[i][1],
			  translations[i][2]};
	vector3_add(trans, vec_frac + 3*i, vec_frac + 3*i);
    }

    // Convert back to cartesian space
    for(int i = 0; i < dim; i++)
	vector3_mat3b3_multiply(lattice_t, vec_frac + 3*i, vec + 3*i);

}

void symmetrize_gradient(double *state, int *invert, const int nmol, const int spg)
{
    // Collect position vectors
    double pos[nmol * 3];
    int dim = 0;

    double lattice[3][3] = {{state[0],        0,       0},
			    {state[1], state[2],       0},
			    {state[3], state[4], state[5]}};

    // each mol is represented by 7 numbers; first 3 are positions
    for(int i = 6; i < 7*nmol + 6; i += 7)
    {
	pos[3*dim + 0] = state[i + 0];
	pos[3*dim + 1] = state[i + 1];
	pos[3*dim + 2] = state[i + 2];
	dim++;
    }
    symmetrize_gradient_pos(pos, dim, lattice, spg);

    // Collect orientational quaternions and convert to rotation matrices
    double quat[4];
    double mat[nmol*3*3];
    dim = 0;
    for(int i = 9; i < 7*nmol + 6; i += 7)
    {
	quat[0] = state[i + 0];
	quat[1] = state[i + 1];
	quat[2] = state[i + 2];
	quat[3] = state[i + 3];
	quaternion2matrix(quat, mat + 3*3*dim, invert[dim]);
	dim++;
    }

    symmetrize_matrix(mat, dim, spg);

    // Convert back to quaternions
    dim = 0;
    double quat_all[nmol*4];
    for(int i = 9; i < 7*nmol + 6; i += 7)
    {
	matrix2quaternion(mat + 3*3*dim, quat_all + 4*dim, invert + dim);
	dim++;
    }

    // Fill in state vector using the results
       
     dim = 0;
    for(int i = 6; i < 7*nmol + 6; i += 7)
    {
	state[i + 0] = pos[3*dim + 0];
	state[i + 1] = pos[3*dim + 1];
	state[i + 2] = pos[3*dim + 2];
	dim ++;
    }
    
    dim = 0;
    for(int i = 9; i < 7*nmol + 6; i += 7)
    {
        state[i + 0] = quat_all[4*dim + 0];
        state[i + 1] = quat_all[4*dim + 1];
        state[i + 2] = quat_all[4*dim + 2];
        state[i + 3] = quat_all[4*dim + 3];
	dim++;
    }
    
    
}



/*
  Symmetrizes a vector wrt symmetry operations of a spg.
  Can be used for symmetrizing positions and position gradients.
  vec, size -> input vector and its dimensions. Shape of vec is (dim, 3).
               dim should be equal to number of symm operation or Z for now
  lattice   -> lattice vectors in row major form
  spg       -> spacegroup which defines the symm operations.
*/
static void symmetrize_vector(double *vec, int dim, double lattice[3][3], int spg)
{
    // Get symm operations
    double translations[192][3];
    int rotations[192][3][3];
    int hall_number = hall_number_from_spg(spg);
    int Z = spg_get_symmetry_from_database(rotations,
					   translations,
					   hall_number);
    assert(dim == Z);
    // Get inverse lattice
    double inverse_lattice[3][3], inverse_lattice_t[3][3], lattice_t[3][3];
    inverse_mat3b3(inverse_lattice, lattice);

    mat3b3_transpose(inverse_lattice_t, inverse_lattice);
    mat3b3_transpose(lattice_t, lattice);

    // Convert to fractional space
    double vec_frac[3*dim];
    for(int i = 0; i < dim; i++)
    {
	vector3_mat3b3_multiply(inverse_lattice_t, vec + 3*i, vec_frac +3*i);
	vector3_frac(vec_frac + 3*i);
    }

    // Compute shifts in fractional space to move molecules
    // to positions where they fall after applying symmetry operation.
    /*
    double shifts[3*dim];
    for(int i = 0; i < dim; i++)
    {
	double temp[3] = {0}, temp_frac[3] = {0};
	int (*rot_i)[3] = rotations[i];
	double rot[3][3] = {{rot_i[0][0], rot_i[0][1],rot_i[0][2]},
			   {rot_i[1][0], rot_i[1][1],rot_i[1][2]},
			   {rot_i[2][0], rot_i[2][1],rot_i[2][2]}};
	vector3_mat3b3_multiply(rot, vec_frac, temp);
	memcpy(temp_frac, temp, 3*sizeof(double));
	vector3_frac(temp_frac);
	vector3_subtract(temp, temp_frac, shifts + 3*i);

	// Shift back to the "natural" position
	vector3_add(shifts + 3*i, vec_frac + 3*i, vec_frac + 3*i);
	
    }
    */
    
    // Array to compute symm average
    double vec_symm[3] = {0};
    
    for(int op = 0; op < Z; op++)
    {
	int (*rot_i)[3] = rotations[op];
	double rot[3][3] = {{rot_i[0][0], rot_i[0][1],rot_i[0][2]},
			   {rot_i[1][0], rot_i[1][1],rot_i[1][2]},
			   {rot_i[2][0], rot_i[2][1],rot_i[2][2]}};
	double trans[3] = {translations[op][0],
			  translations[op][1],
			  translations[op][2]};
	
	double inv_rot[3][3];
	inverse_mat3b3(inv_rot, rot);

	double temp[3];
	vector3_subtract(vec_frac + 3*op, trans, vec_frac + 3*op);
	vector3_mat3b3_multiply(inv_rot, vec_frac + 3*op, temp);
	vector3_frac(temp);
	vector3_add(temp, vec_symm, vec_symm);
    }

    // Take mean over total operations
    for(int i = 0; i < 3; i++)
	vec_symm[i] /= Z;
    
    // Apply symmetry operation to regenerate the entire vec
    for(int i = 0; i < Z; i++)
    {
	int (*rot)[3] = rotations[i];
	vector3_intmat3b3_multiply(rot, vec_symm, vec_frac + 3*i);
	double trans[3] = {translations[i][0],
			  translations[i][1],
			  translations[i][2]};
	vector3_add(trans, vec_frac + 3*i, vec_frac + 3*i);
	vector3_frac(vec_frac + 3*i);
    }

    // Convert back to cartesian space
    for(int i = 0; i < dim; i++)
	vector3_mat3b3_multiply(lattice_t, vec_frac + 3*i, vec + 3*i);
}

void test_symmetrize_vector()
{
    int dim = 2;
    double vec[] = {5, 4, 6, -2, 1, -6};
    double lattice[3][3] = {{5, 0, 0}, {0.5, 3, 0}, {1, 0, 2}};
    for(int i = 0; i < dim; i++)
	printf("%f %f %f\n", vec[0 + 3*i], vec[1 + 3*i], vec[2 + 3*i]);

    symmetrize_vector(vec, dim, lattice, 2);
    for(int i = 0; i < dim; i++)
	printf("%f %f %f\n", vec[0 + 3*i], vec[1 + 3*i], vec[2 + 3*i]);
}

void test_symmetrize_matrix()
{
    double mat[] = { 1.02, 0, 0, 0,  1, 0, 0, 0,  1,
		    -1, 0, 0, 0, -1.05, 0, 0, 0, -1};
    int dim = 2;
    int spg = 2;

   symmetrize_matrix(mat, dim, spg);

    for(int i = 0; i < dim; i++)
    for(int j = 0; j < 3; j++)
	printf("%f, %f, %f \n",
	       mat[9*i + 0 + 3*j],
	       mat[9*i + 1 + 3*j],
	       mat[9*i + 2 + 3*j]
	       );
}

void test_rot_quat_conversion()
{
#define deg2rad(X) X * (3.14/180)
    /*double mat[9] = {1, 0, 0,
		    0, cos(deg2rad(30)), -sin(deg2rad(30)),
		    0, sin(deg2rad(30)), cos(deg2rad(30))};
    */
    double mat[9] = { cos(deg2rad(30)), 0, sin(deg2rad(30)),
			0,              1,              0,
		     -sin(deg2rad(30)), 0, cos(deg2rad(30))};
    double quat[4];
    int inv;
    for(int j = 0; j < 3; j++)
	printf("%f, %f, %f \n",
	       mat[0 + 3*j],
	       mat[1 + 3*j],
	       mat[2 + 3*j]
	       );
    
    matrix2quaternion(mat, quat, &inv);
    printf("quat = %f %f %f %f %d\n", quat[0], quat[1], quat[2], quat[3], inv);

    quaternion2matrix(quat, mat, inv);
    for(int j = 0; j < 3; j++)
	printf("%f, %f, %f \n",
	       mat[0 + 3*j],
	       mat[1 + 3*j],
	       mat[2 + 3*j]
	       );
 }


/*
int main()
{
    test_symmetrize_vector();
    test_symmetrize_matrix();
    //test_rot_quat_conversion();
    test_symmetrize_state();
}
*/
/*
  Symmetrizes a matrix wrt symmetry operations of a spg.
  Useful for rotation matrices.
  mat, size -> input square matrix and its dimension. Shape = (dim, 3, 3)
  lattice   -> lattice vectors in row major form.
  spg       -> spacegroup which defines the symm operations.

  NOTE : Average rotation cannot be calculated using simple mean.
  Assumption is that elements of the input matrices are very close to
  each other.
*/

void symmetrize_matrix(double *mat, int dim, int spg)
{
    // Get symm operations
    double translations[192][3];
    int rotations[192][3][3];
    int hall_number = hall_number_from_spg(spg);
    int Z = spg_get_symmetry_from_database(rotations,
					   translations,
					   hall_number);
    assert(dim == Z);

    // Matrix to compute symm average
    double mat_symm[3][3] = {0};
    double temp[3][3] = {0};
    for(int op = 0; op < Z; op++)
    {
	int (*rot_i)[3] = rotations[op];
	double rot[3][3] = {{rot_i[0][0], rot_i[0][1], rot_i[0][2]},
			   {rot_i[1][0], rot_i[1][1], rot_i[1][2]},
			   {rot_i[2][0], rot_i[2][1], rot_i[2][2]}};

	double inv_rot[3][3]; 
	inverse_mat3b3(inv_rot, rot);

	double (*mat_i)[3] = (double (*)[3]) (mat + 3*3*op);
	mat3b3_mat3b3_multiply(inv_rot, mat_i, temp);
	mat3b3_add(temp, mat_symm, mat_symm);
    }

    // Take average
    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
	mat_symm[i][j] /= Z;

    // For safety, ensure det(R) = 1
    double det = det_mat3b3(mat_symm);
    double renorm = fabsf(1.0 / cbrtf(det));
    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
	mat_symm[i][j] *= renorm;
    
    // Apply symm operation to get full mat
    for(int op = 0; op < Z; op++)
    {
	int (*rot_i)[3] = rotations[op];
	double rot[3][3] = {{rot_i[0][0], rot_i[0][1], rot_i[0][2]},
			   {rot_i[1][0], rot_i[1][1], rot_i[1][2]},
			   {rot_i[2][0], rot_i[2][1], rot_i[2][2]}};

	double (*mat_i)[3] = (double (*)[3]) (mat + 3*3*op);
	mat3b3_mat3b3_multiply(rot, mat_symm, mat_i);
    }
}


// https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
// inv takes -1 or +1 and determines if the rotation is proper or improper
static void quaternion2matrix(double *quat, double *mat, const int inv)
{
    double s = 1.0/(quat[0]*quat[0] +
		   quat[1]*quat[1] +
		   quat[2]*quat[2] +
		   quat[3]*quat[3]);

    mat[0] = 1 - 2 * s * (quat[2]*quat[2] + quat[3]*quat[3]);
    mat[4] = 1 - 2 * s * (quat[1]*quat[1] + quat[3]*quat[3]);
    mat[8] = 1 - 2 * s * (quat[2]*quat[2] + quat[1]*quat[1]);

    mat[3] = 2 * s * (quat[1]*quat[2] - quat[3]*quat[0]);
    mat[6] = 2 * s * (quat[1]*quat[3] + quat[2]*quat[0]);

    mat[1] = 2 * s * (quat[1]*quat[2] + quat[3]*quat[0]);
    mat[7] = 2 * s * (quat[2]*quat[3] - quat[1]*quat[0]);
    
    mat[2] = 2 * s * (quat[1]*quat[3] - quat[2]*quat[0]);
    mat[5] = 2 * s * (quat[2]*quat[3] + quat[1]*quat[0]);

    for(int i = 0; i < 9; i++)
	mat[i] *= inv;
}

// inv -> +1 or -1 depending on proper/improper rotation
static void matrix2quaternion(double *rot, double *quat, int *inv)
{

    double det = rot[0+0*3]*rot[1+1*3]*rot[2+2*3]
	-rot[0+0*3]*rot[1+2*3]*rot[2+1*3]
	+rot[0+1*3]*rot[1+2*3]*rot[2+0*3]
	-rot[0+1*3]*rot[1+0*3]*rot[2+2*3]
	+rot[0+2*3]*rot[1+0*3]*rot[2+1*3]
	-rot[0+2*3]*rot[1+1*3]*rot[2+0*3];

    if(det < 0)
	*inv = -1;
    else
	*inv = +1;

    for(int i = 0; i < 9; i++)
	rot[i] *= *inv;
    
    double t;
    if (rot[2+2*3] < 0.0)
    {
        if (rot[0+0*3] > rot[1+1*3])
        {
            t = 1.0 + rot[0+0*3] - rot[1+1*3] - rot[2+2*3];
            quat[0] = rot[2+1*3]-rot[1+2*3];
            quat[1] = t;
            quat[2] = rot[1+0*3]+rot[0+1*3];
            quat[3] = rot[2+0*3]+rot[0+2*3];
        }
        else
        {
            t = 1.0 - rot[0+0*3] + rot[1+1*3] - rot[2+2*3];
            quat[0] = rot[0+2*3]-rot[2+0*3];
            quat[1] = rot[1+0*3]+rot[0+1*3];
            quat[2] = t;
            quat[3] = rot[2+1*3]+rot[1+2*3];
        }
    }
    else
    {
        if (rot[0+0*3] < -rot[1+1*3])
        {
            t = 1.0 - rot[0+0*3] - rot[1+1*3] + rot[2+2*3];
            quat[0] = rot[1+0*3]-rot[0+1*3];
            quat[1] = rot[0+2*3]+rot[2+0*3];
            quat[2] = rot[2+1*3]+rot[1+2*3];
            quat[3] = t;
        }
        else
        {
            t = 1.0 + rot[0+0*3] + rot[1+1*3] + rot[2+2*3];
            quat[0] = t;
            quat[1] = rot[2+1*3]-rot[1+2*3];
            quat[2] = rot[0+2*3]-rot[2+0*3];
            quat[3] = rot[1+0*3]-rot[0+1*3];
        }
    }
    t = 0.5/sqrt(t);
    quat[0] *= t;
    quat[1] *= t;
    quat[2] *= t;
    quat[3] *= t;
    
}

void test_symmetrize_state()
{
    double state[6 + 2*7] = {1, 0, 1, 0, 0, 1,
			     5, 4, 6,
			     1, 0, 0, 0,   
			     -2, 1, -6,
			     1, 0, 0, 0};   
    int nmol = 2;
    int invert[2] = {1, -1};
    int spg = 2;
    symmetrize_state(state, invert, nmol, spg);
    for(int i = 0; i < 7*nmol + 6; i++)
	printf("%f ", state[i]);
}

void symmetrize_state(double *state, int *invert, const int nmol, const int spg)
{
    // Collect position vectors
    double pos[nmol * 3];
    int dim = 0;

    double lattice[3][3] = {{state[0],        0,       0},
			    {state[1], state[2],       0},
			    {state[3], state[4], state[5]}};

    // each mol is represented by 7 numbers; first 3 are positions
    for(int i = 6; i < 7*nmol + 6; i += 7)
    {
	pos[3*dim + 0] = state[i + 0];
	pos[3*dim + 1] = state[i + 1];
	pos[3*dim + 2] = state[i + 2];
	dim++;
    }
    symmetrize_vector(pos, dim, lattice, spg);

    // Collect orientational quaternions and convert to rotation matrices
    double quat[4];
    double mat[nmol*3*3];
    dim = 0;
    for(int i = 9; i < 7*nmol + 6; i += 7)
    {
	quat[0] = state[i + 0];
	quat[1] = state[i + 1];
	quat[2] = state[i + 2];
	quat[3] = state[i + 3];
	quaternion2matrix(quat, mat + 3*3*dim, invert[dim]);
	dim++;
    }

    symmetrize_matrix(mat, dim, spg);

    // Convert back to quaternions
    dim = 0;
    double quat_all[nmol*4];
    for(int i = 9; i < 7*nmol + 6; i += 7)
    {
	matrix2quaternion(mat + 3*3*dim, quat_all + 4*dim, invert + dim);
	dim++;
    }

    // Fill in state vector using the results
       
     dim = 0;
    for(int i = 6; i < 7*nmol + 6; i += 7)
    {
	state[i + 0] = pos[3*dim + 0];
	state[i + 1] = pos[3*dim + 1];
	state[i + 2] = pos[3*dim + 2];
	dim ++;
    }
    
    dim = 0;
    for(int i = 9; i < 7*nmol + 6; i += 7)
    {
        state[i + 0] = quat_all[4*dim + 0];
        state[i + 1] = quat_all[4*dim + 1];
        state[i + 2] = quat_all[4*dim + 2];
        state[i + 3] = quat_all[4*dim + 3];
	dim++;
    }
    
    
}

int detect_spg_using_spglib(crystal* xtal)
{
    double tol = 1;
    //print_crystal(xtal)
    convert_xtal_to_fractional(xtal);
    //variable declarations
    int num_atoms_in_cell = xtal->Z * xtal-> num_atoms_in_molecule ;
    int types[num_atoms_in_cell];
    double positions[num_atoms_in_cell][3];
    char atom[num_atoms_in_cell*2];
    char symbol[21];
    int num_spg;
    double lattice_vector[3][3];
    //printf("total atoms = %d \n", num_atoms_in_cell);
    for(int i = 0; i < num_atoms_in_cell; i++)
    {
        positions[i][0] = xtal->Xcord[i];
        positions[i][1] = xtal->Ycord[i];
        positions[i][2] = xtal->Zcord[i];
        atom[2*i] = xtal->atoms[2*i];
        atom[2*i+1] = xtal->atoms[2*i+1];

        if      (atom[2*i] == 'C' )
        types[i] = 6;
        else if (atom[2*i] == 'H' )
        types[i] = 1;
        else if (atom[2*i] == 'N' && atom[2*i+1] == ' ')
        types[i] = 7;
        else if (atom[2*i] == 'O' && atom[2*i+1] == ' ')
        types[i] = 8;
        else if (atom[2*i] == 'F' && atom[2*i+1] == ' ')
        types[i] = 9;
        else if (atom[2*i] == 'P' && atom[2*i+1] == ' ')
        types[i] = 15;
        else if (atom[2*i] == 'S' && atom[2*i+1] == ' ')
        types[i] = 16;
        else if (atom[2*i] == 'B' && atom[2*i+1] == 'r')
        types[i] = 35;
        else if (atom[2*i] == 'I' && atom[2*i+1] == ' ')
        types[i] = 53;
        else if (atom[2*i] == 'B' && atom[2*i+1] == ' ')
        types[i] = 5;
        else if (atom[2*i] == 'H' && atom[2*i+1] == 'e')
        types[i] = 2;
        else if (atom[2*i] == 'N' && atom[2*i+1] == 'e')
        types[i] = 10;
        else if (atom[2*i] == 'K' && atom[2*i+1] == 'r')
        types[i] = 36;
        else if (atom[2*i] == 'S' && atom[2*i+1] == 'i')
        types[i] = 14;
        else
        {printf("***ERROR: spglib detector: atom not found -> %c%c\n", atom[2*i], atom[2*i+1]);exit(0);}

    }
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++ )
        {
            //note transpose
            lattice_vector[i][j] = xtal->lattice_vectors[j][i];
        }

    num_spg = spg_get_international(symbol,
                                    lattice_vector,
                                    positions,
                                    types,
                                    num_atoms_in_cell,
                                    tol);
    //SpglibError error;
    //error = spg_get_error_code();
    //printf("#SPGLIB says %s\n", spg_get_error_message(error));
    //printf("#SPGLIB detected space group = %d\n",
    //                                         num_spg);
    convert_xtal_to_cartesian(xtal);
    return num_spg;

}

void convert_xtal_to_cartesian(crystal *xtal)
{
    double lattice_vectors_transpose[3][3] =
	{
	    {xtal->lattice_vectors[0][0], xtal->lattice_vectors[0][1], xtal->lattice_vectors[0][2]},
	    {xtal->lattice_vectors[1][0], xtal->lattice_vectors[1][1], xtal->lattice_vectors[1][2]},
	    {xtal->lattice_vectors[2][0], xtal->lattice_vectors[2][1], xtal->lattice_vectors[2][2]}
	};

    mat3b3_transpose(lattice_vectors_transpose,
                     lattice_vectors_transpose);
    int total_atoms = xtal->Z *xtal->num_atoms_in_molecule;

    for(int i = 0; i < total_atoms; i++)
    {
        double temp[3] = {xtal->Xcord[i],
                         xtal->Ycord[i],
                         xtal->Zcord[i] };
        vector3_mat3b3_multiply(lattice_vectors_transpose, temp, temp );
        xtal->Xcord[i] = temp[0];
        xtal->Ycord[i] = temp[1];
        xtal->Zcord[i] = temp[2];
    }
}

void convert_xtal_to_fractional(crystal *xtal)
{
    double lattice_vectors[3][3] ={
	{xtal->lattice_vectors[0][0], xtal->lattice_vectors[0][1], xtal->lattice_vectors[0][2]},
	{xtal->lattice_vectors[1][0], xtal->lattice_vectors[1][1], xtal->lattice_vectors[1][2]},
	{xtal->lattice_vectors[2][0], xtal->lattice_vectors[2][1], xtal->lattice_vectors[2][2]}
    };

    double inv_lattice_vectors[3][3];
    inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
    mat3b3_transpose(inv_lattice_vectors, inv_lattice_vectors);

    int total_atoms = xtal->Z *xtal->num_atoms_in_molecule;
    
    for(int i = 0; i < total_atoms; i++)
    {
        double temp[3] = {xtal->Xcord[i], xtal->Ycord[i], xtal->Zcord[i] };
        vector3_mat3b3_multiply(inv_lattice_vectors, temp, temp );
        xtal->Xcord[i] = temp[0];
        xtal->Ycord[i] = temp[1];
        xtal->Zcord[i] = temp[2];
    }
}

void bring_all_molecules_to_first_cell(crystal *xtal)
{
    int N = xtal->num_atoms_in_molecule;
    int Z = xtal->Z;

    convert_xtal_to_fractional(xtal);
    //convert to fractional
    for(int i = 0; i < N*Z; i+=N)
    {
        double com[3] = {0,0,0};
        compute_molecule_COM(*xtal,com,i);
        int integer_part1 = (int)com[0];
        int integer_part2 = (int)com[1];
        int integer_part3 = (int)com[2];
        if(com[0] < 0)
        {integer_part1 -= 1;}
        if(com[1] < 0)
        {integer_part2 -= 1;}
        if(com[2] < 0)
        {integer_part3 -= 1;}

        for(int j = i; j < i + N; j++)
        {
            xtal->Xcord[j] -= integer_part1;
            xtal->Ycord[j] -= integer_part2;
            xtal->Zcord[j] -= integer_part3;
        }

    }

    convert_xtal_to_cartesian(xtal);
}

void compute_molecule_COM(crystal xtal, double com[3], int i)
{   // computes the COM of a molecule (whose first atom is i) in
    //a crystal.
    
    com[0] = 0;
    com[1] = 0;
    com[2] = 0;

    int N = xtal.num_atoms_in_molecule;
    for(int j = i; j < i + N ; j++)
    {
        com[0] += xtal.Xcord[j];
        com[1] += xtal.Ycord[j];
        com[2] += xtal.Zcord[j];
    }
    com[0] /= N;
    com[1] /= N;
    com[2] /= N;
}

void print_crystal(crystal* xtal)
{
     int N = xtal->num_atoms_in_molecule;
     int m = xtal->Z;

     printf("#Z = %d \n", xtal->Z);
     printf("#napm = %d \n", xtal->num_atoms_in_molecule);

    for(int i = 0; i < 3; i++)
    {
        printf("lattice_vector %12f %12f %12f \n",
            xtal->lattice_vectors[i][0], xtal->lattice_vectors[i][1],
            xtal->lattice_vectors[i][2]);
    }

    for(int i = 0; i < N*m; i++)
    {
        printf("atom %12f %12f %12f  %c \n", xtal->Xcord[i],
            xtal->Ycord[i],  xtal->Zcord[i],  xtal->atoms[2*i]);
    }
}


