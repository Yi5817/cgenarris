// Geometry optimization of molecular crystals in a regularized rigid-body approximation

// METHODOLOGY COMMENTS ON THIS OPTIMIZER:
// - It would be more efficient & robust to work with only the symmetry-inequivalent structural degrees of freedom,
//   but we don't have a convenient way of expressing those degrees of freedom, especially for Wyckoff positions.
//
// - It would be a little more efficient & stable to only consider 3 explicit degrees of freedom when perturbing
//   the orientational quaternions. For example, rather than an additive perturbation, a multiplicative perturbation
//   of the form (1, a, b, c) could be considered for infinitesimal a, b, & c.
//
// - There is nothing particularly special about the regularized rigid-body approximation that was chosen here,
//   it could be replaced by a more physical & accurate intermolecular interaction, which might improve the overlap
//   between the set of structures that are local minima here versus the true local minima of an accurate energy surface.
//
// - Space-group symmetries are likely to break if the initial structure is too loosely packed relative to the cutoff
//   distance of the regularized interatomic interaction. This can be mitigated by increasing the cutoff distance.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../../crystal.h"
#include "../../crystal_utils.h"
#include "rigid_press.h"
#include "symmetrization.h"

//prototypes for external dependencies (BLAS & LAPACK)
void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
void dgeqrf_(int*, int*, double*, int*, double*, double*, int*, int*);
void dormqr_(char*, char*, int*, int*, int*, double*, int*, double*, double*, int*, double*, int*, int*);
void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int* ,int*);

// cutoff distance of the interaction kernel
#define INTERACTION_CUTOFF 10.0

// parameters defining the regularized interatomic contact interaction
#define INTERACTION_WEIGHT 0.1

// maximum number of unit cells to sum over before deciding a crystal is too tightly packed
#define MAX_CELL_SUM 1000

// number of steps to take for each Golden-section line search
#define GOLDEN_STEPS 20

// optimization tolerance on energy change in 1 iteration relative to minimum energy
#define OPTIMIZATION_TOLERANCE 1e-8

// approximate maximum 2-norm change in the state vector for line searches (scaled by nmol)
#define MAX_LINE_STEP 100.0

// very loose energy reduction threshold for flagging false convergence (scaled by nmol)
#define ENERGY_REDUCTION_TOLERANCE 1e-3

// step size for numerical tests of analytical derivatives
#define STEP 1e-4

// the state vector of a crystal's geometry has size 6+7*nmol containing the following information:
// state[0] : 1st lattice vector (x)
// state[1-2] : 2nd lattice vector (x-y)
// state[3-5] : 3rd lattice vector (x-y-z)
// state[6-8] : center of 1st molecule (x-y-z)
// state[9-12] : orientational quaternion of 1st molecule
// ...

// information about a molecular crystal to define the structural relaxation problem
struct molecular_crystal
{
    // size information
    int ntype; // number of types of molecules in the crystal
    int nmol; // number of molecules per unit cell
    int spg;
    // information about each type of molecule
    int *natom; // number of atoms in a molecule type [ntype]
    double **geometry; // centered geometry of a molecule type (interleaved x-y-z format) [ntype][3*natom[i]]
    double ***collide; // collision distances for pairs of atoms between 2 molecules [ntype][ntype][natom[i]*natom[j]]

    // information about the crystal
    double *mol_length; // Length of molecule. distance from geometric center to the farthest atom[nmol]
    int *type; // type of each molecule in the unit cell [nmol]
    int *invert; // inversion of each molecule in the unit cell (+1 for standard, -1 for inverted) [nmol]
};

#ifdef ROPT_DEBUG
void state_2_xtal(crystal *xtl, double *state, struct molecular_crystal *xtl2);
void print_crystal(crystal* xtal);
int detect_spg_using_spglib(crystal* xtal);
#endif
static int get_cell_type_from_spg(int spg);


// deallocate memory in the molecular_crystal structure
void free_molecular_crystal(struct molecular_crystal *xtl)
{
    free(xtl->type);
    free(xtl->invert);
    free(xtl->mol_length);

    for(int i=0 ; i<xtl->ntype ; i++)
    { free(xtl->geometry[i]); }
    free(xtl->geometry);

    for(int i=0 ; i<xtl->ntype ; i++)
    {
        for(int j=0 ; j<xtl->ntype ; j++)
        { free(xtl->collide[i][j]); }
	free(xtl->collide[i]);
    }
    free(xtl->collide);

    free(xtl->natom);
}

#ifdef ROPT_DEBUG
void print_state(double *state, int size)
{
    for(int i = 6; i < size; i++)
	printf("%lf ", state[i]);
    printf("\n\n");
}
#endif

// interatomic interaction kernel & its 1st & 2nd distance derivatives
static inline double kernel(double distance, // interatomic distance
              // kernel parameters: (add/replace parameters for more physical interactions)
              double distance0, // collision distance
              double wt) // kernel weight
{
    if(distance > INTERACTION_CUTOFF)
    { return 0.0; }
    if(distance < distance0)
    { return INFINITY; }
    return wt*(INTERACTION_CUTOFF - distance)/(distance - distance0);
}
static inline double dkernel(double distance, // interatomic distance
              // kernel parameters: (add/replace parameters for more physical interactions)
              double distance0, // collision distance
              double wt) // kernel weight
{
    if(distance > INTERACTION_CUTOFF)
    { return 0.0; }
    if(distance < distance0)
    { return NAN; }
    double recip = 1.0/(distance - distance0);
    return -wt*((INTERACTION_CUTOFF - distance)*recip*recip + recip);
}
static inline double d2kernel(double distance, // interatomic distance
              // kernel parameters: (add/replace parameters for more physical interactions)
              double distance0, // collision distance
              double wt) // kernel weight
{
    if(distance > INTERACTION_CUTOFF)
    { return 0.0; }
    if(distance < distance0)
    { return NAN; }
    double recip = 1.0/(distance - distance0);
    return 2.0*wt*((INTERACTION_CUTOFF - distance)*recip*recip*recip + recip*recip);
}

// position of atom in a translated & rotated molecule
static inline void position(double *restrict local, // local coordinate of atom in the molecule [3]
              double *restrict state, // state vector of the molecule (x-y-z & orientational quaternion) [7]
              double *restrict global) // global coordinate of the atom in the crystal [3]
{
    double wt = 2.0/(state[3]*state[3] + state[4]*state[4] + state[5]*state[5] + state[6]*state[6]);
    global[0] = state[0] + local[0] + wt*(-(state[5]*state[5] + state[6]*state[6])*local[0]
                                         + (state[4]*state[5] - state[3]*state[6])*local[1]
                                         + (state[4]*state[6] + state[3]*state[5])*local[2]);
    global[1] = state[1] + local[1] + wt*( (state[4]*state[5] + state[3]*state[6])*local[0]
                                         - (state[4]*state[4] + state[6]*state[6])*local[1]
                                         + (state[5]*state[6] - state[3]*state[4])*local[2]);
    global[2] = state[2] + local[2] + wt*( (state[4]*state[6] - state[3]*state[5])*local[0]
                                         + (state[5]*state[6] + state[3]*state[4])*local[1]
                                         - (state[4]*state[4] + state[5]*state[5])*local[2]);
}
// NOTE: simple derivatives w.r.t. state[0], state[1], & state[2] are ignored here
void position_derivative(double *restrict local, // local coordinate of atom in the molecule [3]
                         double *restrict state, // state vector of the molecule (x-y-z & orientational quaternion) [7]
                         double *restrict global1, // global coordinate 1st derivatives of the atom in the crystal [12]
                         double *restrict global2) // global coordinate 2nd derivatives of the atom in the crystal [48]
{
    double wt = 2.0/(state[3]*state[3] + state[4]*state[4] + state[5]*state[5] + state[6]*state[6]);
    double dwt = -wt*wt, dglobal[3], dglobal2[12];
    dglobal[0] = dwt*(-(state[5]*state[5] + state[6]*state[6])*local[0]
                     + (state[4]*state[5] - state[3]*state[6])*local[1]
                     + (state[4]*state[6] + state[3]*state[5])*local[2]);
    dglobal[1] = dwt*( (state[4]*state[5] + state[3]*state[6])*local[0]
                     - (state[4]*state[4] + state[6]*state[6])*local[1]
                     + (state[5]*state[6] - state[3]*state[4])*local[2]);
    dglobal[2] = dwt*( (state[4]*state[6] - state[3]*state[5])*local[0]
                     + (state[5]*state[6] + state[3]*state[4])*local[1]
                     - (state[4]*state[4] + state[5]*state[5])*local[2]);
    dglobal2[0 + 0*3] = dwt*(-state[6]*local[1] + state[5]*local[2]);
    dglobal2[0 + 1*3] = dwt*( state[5]*local[1] + state[6]*local[2]);
    dglobal2[0 + 2*3] = dwt*(-2.0*state[5]*local[0] + state[4]*local[1] + state[3]*local[2]);
    dglobal2[0 + 3*3] = dwt*(-2.0*state[6]*local[0] - state[3]*local[1] + state[4]*local[2]);

    dglobal2[1 + 0*3] = dwt*(state[6]*local[0] - state[4]*local[2]);
    dglobal2[1 + 1*3] = dwt*(state[5]*local[0] - 2.0*state[4]*local[1] - state[3]*local[2]);
    dglobal2[1 + 2*3] = dwt*(state[4]*local[0] + state[6]*local[2]);
    dglobal2[1 + 3*3] = dwt*(state[3]*local[0] - 2.0*state[6]*local[1] + state[5]*local[2]);

    dglobal2[2 + 0*3] = dwt*(-state[5]*local[0] + state[4]*local[1]);
    dglobal2[2 + 1*3] = dwt*( state[6]*local[0] + state[3]*local[1] - 2.0*state[4]*local[2]);
    dglobal2[2 + 2*3] = dwt*(-state[3]*local[0] + state[6]*local[1] - 2.0*state[5]*local[2]);
    dglobal2[2 + 3*3] = dwt*( state[4]*local[0] + state[5]*local[1]);

    for(int i=0 ; i<4 ; i++)
    for(int j=0 ; j<3 ; j++)
    { global1[j + i*3] = state[i+3]*dglobal[j]; }

    global1[0 + 0*3] += wt*(-state[6]*local[1] + state[5]*local[2]);
    global1[0 + 1*3] += wt*( state[5]*local[1] + state[6]*local[2]);
    global1[0 + 2*3] += wt*(-2.0*state[5]*local[0] + state[4]*local[1] + state[3]*local[2]);
    global1[0 + 3*3] += wt*(-2.0*state[6]*local[0] - state[3]*local[1] + state[4]*local[2]);

    global1[1 + 0*3] += wt*(state[6]*local[0] - state[4]*local[2]);
    global1[1 + 1*3] += wt*(state[5]*local[0] - 2.0*state[4]*local[1] - state[3]*local[2]);
    global1[1 + 2*3] += wt*(state[4]*local[0] + state[6]*local[2]);
    global1[1 + 3*3] += wt*(state[3]*local[0] - 2.0*state[6]*local[1] + state[5]*local[2]);

    global1[2 + 0*3] += wt*(-state[5]*local[0] + state[4]*local[1]);
    global1[2 + 1*3] += wt*( state[6]*local[0] + state[3]*local[1] - 2.0*state[4]*local[2]);
    global1[2 + 2*3] += wt*(-state[3]*local[0] + state[6]*local[1] - 2.0*state[5]*local[2]);
    global1[2 + 3*3] += wt*( state[4]*local[0] + state[5]*local[1]);

    for(int i=0 ; i<4 ; i++)
    for(int j=0 ; j<4 ; j++)
    for(int k=0 ; k<3 ; k++)
    {
        global2[k + j*3 + i*12] = -2.0*wt*state[i+3]*state[j+3]*dglobal[k]
                                + state[i+3]*dglobal2[k + j*3] + state[j+3]*dglobal2[k + i*3];
    }

    for(int i=0 ; i<4 ; i++)
    for(int j=0 ; j<3 ; j++)
    { global2[j + i*3 + i*12] += dglobal[j]; }

    global2[0 + 0*3 + 2*12] += wt*local[2];
    global2[0 + 0*3 + 3*12] += -wt*local[1];
    global2[0 + 1*3 + 2*12] += wt*local[1];
    global2[0 + 1*3 + 3*12] += wt*local[2];
    global2[0 + 2*3 + 0*12] += wt*local[2];
    global2[0 + 2*3 + 1*12] += wt*local[1];
    global2[0 + 2*3 + 2*12] += -2.0*wt*local[0];
    global2[0 + 3*3 + 0*12] += -wt*local[1];
    global2[0 + 3*3 + 1*12] += wt*local[2];
    global2[0 + 3*3 + 3*12] += -2.0*wt*local[0];

    global2[1 + 0*3 + 1*12] += -wt*local[2];
    global2[1 + 0*3 + 3*12] += wt*local[0];
    global2[1 + 1*3 + 0*12] += -wt*local[2];
    global2[1 + 1*3 + 1*12] += -2.0*wt*local[1];
    global2[1 + 1*3 + 2*12] += wt*local[0];
    global2[1 + 2*3 + 1*12] += wt*local[0];
    global2[1 + 2*3 + 3*12] += wt*local[2];
    global2[1 + 3*3 + 0*12] += wt*local[0];
    global2[1 + 3*3 + 2*12] += wt*local[2];
    global2[1 + 3*3 + 3*12] += -2.0*wt*local[1];

    global2[2 + 0*3 + 1*12] += wt*local[1];
    global2[2 + 0*3 + 2*12] += -wt*local[0];
    global2[2 + 1*3 + 0*12] += wt*local[1];
    global2[2 + 1*3 + 1*12] += -2.0*wt*local[2];
    global2[2 + 1*3 + 3*12] += wt*local[0];
    global2[2 + 2*3 + 0*12] += -wt*local[0];
    global2[2 + 2*3 + 2*12] += -2.0*wt*local[2];
    global2[2 + 2*3 + 3*12] += wt*local[1];
    global2[2 + 3*3 + 1*12] += wt*local[0];
    global2[2 + 3*3 + 2*12] += wt*local[1];
}
void position_derivative_test(double *local, // local coordinate of atom in the molecule [3]
                              double *state, // state vector of the molecule [7]
                              double *global1, // global coordinate 1st derivatives of the atom in the crystal [12]
                              double *global2) // global coordinate 2nd derivatives of the atom in the crystal [48]
{
    for(int i=0 ; i<4 ; i++)
    {
        double state0 = state[i+3], global_plus[3], global_minus[3], global1_plus[12], global1_minus[12], dummy[48];
        state[i+3] += STEP;
        position(local, state, global_plus);
        position_derivative(local, state, global1_plus, dummy);
        state[i+3] -= 2.0*STEP;
        position(local, state, global_minus);
        position_derivative(local, state, global1_minus, dummy);
        state[i+3] = state0;

        for(int j=0 ; j<3 ; j++)
        { global1[j + i*3] = (global_plus[j] - global_minus[j])/(2.0*STEP); }
        for(int j=0 ; j<12 ; j++)
        { global2[j + i*12] = (global1_plus[j] - global1_minus[j])/(2.0*STEP); }
    }
}

// interaction energy between a pair of molecules
double pair_energy(int natom1, // number of atoms in the 1st molecule
                   int natom2, // number of atoms in the 2nd molecule
                   int invert1, // inversion flag for 1st molecule {+1, -1}
                   int invert2, // inversion flag for 2nd molecule {+1, -1}
                   double *geo1, // geometry of the 1st molecule [3*natom1]
                   double *geo2, // geometry of the 2nd molecule [3*natom2]
                   double *collide, // collision matrix [natom1*natom2]
                   double *state1, // state vector of the 1st molecule [7]
                   double *state2) // state vector of the 2nd molecule [7]
{
    double energy = 0.0, wt = INTERACTION_WEIGHT/(double)(natom1*natom2);
    for(int i=0 ; i<natom1 ; i++)
    {
        double local1[3], coord1[3];
        local1[0] = invert1*geo1[3*i]; local1[1] = invert1*geo1[1+3*i]; local1[2] = invert1*geo1[2+3*i];
        position(local1, state1, coord1);

        for(int j=0 ; j<natom2 ; j++)
        {
            double local2[3], coord2[3];
            local2[0] = invert2*geo2[3*j]; local2[1] = invert2*geo2[1+3*j]; local2[2] = invert2*geo2[2+3*j];
            position(local2, state2, coord2);

            // regularized contact interaction
            double dist = sqrt((coord1[0] - coord2[0])*(coord1[0] - coord2[0])
                              +(coord1[1] - coord2[1])*(coord1[1] - coord2[1])
                              +(coord1[2] - coord2[2])*(coord1[2] - coord2[2]));

            energy += kernel(dist, collide[i+j*natom1], wt);
	    if(energy == INFINITY)
	    { return energy; }
        }
    }
    return energy;
}
void pair_energy_derivative(int natom1, // number of atoms in the 1st molecule
                            int natom2, // number of atoms in the 2nd molecule
                            int invert1, // inversion flag for 1st molecule {+1, -1}
                            int invert2, // inversion flag for 2nd molecule {+1, -1}
                            double *geo1, // geometry of the 1st molecule [3*natom1]
                            double *geo2, // geometry of the 2nd molecule [3*natom2]
                            double *collide, // collision matrix [natom1*natom2]
                            double *state1, // state vector of the 1st molecule [7]
                            double *state2, // state vector of the 2nd molecule [7]
                            double *grad1, // 1st derivatives w.r.t. 1st state vector [7]
                            double *grad2, // 1st derivatives w.r.t. 2nd state vector [7]
                            double *hess11, // 2nd derivatives w.r.t. 1st state vector [49]
                            double *hess22, // 2nd derivatives w.r.t. 2nd state vector [49]
                            double *hess12) // mixed 2nd derivatives [49]
{
    for(int i=0 ; i<7 ; i++)
    { grad1[i] = grad2[i] = 0.0; }

    for(int i=0 ; i<49 ; i++)
    { hess11[i] = hess12[i] = hess22[i] = 0.0; }

    double wt = INTERACTION_WEIGHT/(double)(natom1*natom2);
    for(int i=0 ; i<natom1 ; i++)
    {
        double local1[3], coord1[3], coord1_deriv1[12], coord1_deriv2[48];
        local1[0] = invert1*geo1[3*i]; local1[1] = invert1*geo1[1+3*i]; local1[2] = invert1*geo1[2+3*i];
        position(local1, state1, coord1);
        position_derivative(local1, state1, coord1_deriv1, coord1_deriv2);

        for(int j=0 ; j<natom2 ; j++)
        {
            double local2[3], coord2[3], coord2_deriv1[12], coord2_deriv2[48];
            local2[0] = invert2*geo2[3*j]; local2[1] = invert2*geo2[1+3*j]; local2[2] = invert2*geo2[2+3*j];
            position(local2, state2, coord2);
            position_derivative(local2, state2, coord2_deriv1, coord2_deriv2);

            double dist = sqrt((coord1[0] - coord2[0])*(coord1[0] - coord2[0])
                              +(coord1[1] - coord2[1])*(coord1[1] - coord2[1])
                              +(coord1[2] - coord2[2])*(coord1[2] - coord2[2]));

	    // If the atoms don't interact, derivatives are trivially zero.
	    if(dist > INTERACTION_CUTOFF)
	    { continue; }
	    
            double denergy = dkernel(dist, collide[i+j*natom1], wt)/dist;
            double d2energy = (d2kernel(dist, collide[i+j*natom1], wt) - denergy)/(dist*dist);

            // simple derivatives w.r.t. state1[0], state1[1], state1[2], state2[0], state2[1], & state2[2]
            for(int k=0 ; k<3 ; k++)
            {
                double delta = coord1[k] - coord2[k];
                grad1[k] += denergy*delta;
                grad2[k] -= denergy*delta;
                hess11[k + k*7] += denergy;
                hess22[k + k*7] += denergy;
                hess12[k + k*7] -= denergy;
                for(int l=0 ; l<3 ; l++)
                {
                    double delta2 = coord1[l] - coord2[l];
                    hess11[l + k*7] += d2energy*delta*delta2;
                    hess22[l + k*7] += d2energy*delta*delta2;
                    hess12[l + k*7] -= d2energy*delta*delta2;
                }
            }

            // gradients & mixed quaternion-translation hessian
            for(int k=0 ; k<4 ; k++)
            for(int l=0 ; l<3 ; l++)
            {
                double delta = coord1[l] - coord2[l];
                grad1[k+3] += denergy*delta*coord1_deriv1[l+k*3];
                grad2[k+3] -= denergy*delta*coord2_deriv1[l+k*3];
                hess11[k+3 + l*7] += denergy*coord1_deriv1[l+k*3];
                hess22[k+3 + l*7] += denergy*coord2_deriv1[l+k*3];
                hess12[k+3 + l*7] -= denergy*coord1_deriv1[l+k*3];
                hess11[l + (k+3)*7] += denergy*coord1_deriv1[l+k*3];
                hess22[l + (k+3)*7] += denergy*coord2_deriv1[l+k*3];
                hess12[l + (k+3)*7] -= denergy*coord2_deriv1[l+k*3];

                for(int m=0 ; m<3 ; m++)
                {
                    double delta2 = coord1[m] - coord2[m];
                    hess11[k+3 + m*7] += d2energy*delta*delta2*coord1_deriv1[l+k*3];
                    hess22[k+3 + m*7] += d2energy*delta*delta2*coord2_deriv1[l+k*3];
                    hess12[k+3 + m*7] -= d2energy*delta*delta2*coord1_deriv1[l+k*3];
                    hess11[m + (k+3)*7] += d2energy*delta*delta2*coord1_deriv1[l+k*3];
                    hess22[m + (k+3)*7] += d2energy*delta*delta2*coord2_deriv1[l+k*3];
                    hess12[m + (k+3)*7] -= d2energy*delta*delta2*coord2_deriv1[l+k*3];
                }
            }

            // quaternion components of hessian
            for(int k=0 ; k<4 ; k++)
            for(int l=0 ; l<4 ; l++)
            for(int m=0 ; m<3 ; m++)
            {
                double delta = coord1[m] - coord2[m];
                hess11[k+3 + (l+3)*7] += denergy*coord1_deriv1[m+k*3]*coord1_deriv1[m+l*3];
                hess22[k+3 + (l+3)*7] += denergy*coord2_deriv1[m+k*3]*coord2_deriv1[m+l*3];
                hess12[k+3 + (l+3)*7] -= denergy*coord1_deriv1[m+k*3]*coord2_deriv1[m+l*3];
                hess11[k+3 + (l+3)*7] += denergy*delta*coord1_deriv2[m+k*3+l*12];
                hess22[k+3 + (l+3)*7] -= denergy*delta*coord2_deriv2[m+k*3+l*12];

                for(int n=0 ; n<3 ; n++)
                {
                    double delta2 = coord1[n] - coord2[n];
                    hess11[k+3 + (l+3)*7] += d2energy*delta*delta2*coord1_deriv1[m+k*3]*coord1_deriv1[n+l*3];
                    hess22[k+3 + (l+3)*7] += d2energy*delta*delta2*coord2_deriv1[m+k*3]*coord2_deriv1[n+l*3];
                    hess12[k+3 + (l+3)*7] -= d2energy*delta*delta2*coord1_deriv1[m+k*3]*coord2_deriv1[n+l*3];
                }
            }
        }
    }
}
void pair_energy_derivative_test(int natom1, // number of atoms in the 1st molecule
                                 int natom2, // number of atoms in the 2nd molecule
                                 int invert1, // inversion flag for 1st molecule {+1, -1}
                                 int invert2, // inversion flag for 2nd molecule {+1, -1}
                                 double *geo1, // geometry of the 1st molecule [3*natom1]
                                 double *geo2, // geometry of the 2nd molecule [3*natom2]
                                 double *collide, // collision matrix [natom1*natom2]
                                 double *state1, // state vector of the 1st molecule [7]
                                 double *state2, // state vector of the 2nd molecule [7]
                                 double *grad1, // 1st derivatives w.r.t. 1st state vector [7]
                                 double *grad2, // 1st derivatives w.r.t. 2nd state vector [7]
                                 double *hess11, // 2nd derivatives w.r.t. 1st state vector [49]
                                 double *hess22, // 2nd derivatives w.r.t. 2nd state vector [49]
                                 double *hess12) // mixed 2nd derivatives [49]
{
    for(int i=0 ; i<7; i++)
    {
        double state0 = state1[i], grad1_plus[7], grad1_minus[7], grad2_plus[7], grad2_minus[7], dummy[49];
        state1[i] += STEP;
        double energy_plus = pair_energy(natom1,natom2,invert1,invert2,geo1,geo2,collide,state1,state2);
        pair_energy_derivative(natom1,natom2,invert1,invert2,geo1,geo2,collide,state1,state2,grad1_plus,grad2_plus,dummy,dummy,dummy);
        state1[i] -= 2.0*STEP;
        double energy_minus = pair_energy(natom1,natom2,invert1,invert2,geo1,geo2,collide,state1,state2);
        pair_energy_derivative(natom1,natom2,invert1,invert2,geo1,geo2,collide,state1,state2,grad1_minus,grad2_minus,dummy,dummy,dummy);
        state1[i] = state0;

        grad1[i] = (energy_plus - energy_minus)/(2.0*STEP);
        for(int j=0 ; j<7 ; j++)
        {
            hess11[i+j*7] = (grad1_plus[j] - grad1_minus[j])/(2.0*STEP);
            hess12[i+j*7] = (grad2_plus[j] - grad2_minus[j])/(2.0*STEP);
        }

        state0 = state2[i];
        state2[i] += STEP;
        energy_plus = pair_energy(natom1,natom2,invert1,invert2,geo1,geo2,collide,state1,state2);
        pair_energy_derivative(natom1,natom2,invert1,invert2,geo1,geo2,collide,state1,state2,grad1_plus,grad2_plus,dummy,dummy,dummy);
        state2[i] -= 2.0*STEP;
        energy_minus = pair_energy(natom1,natom2,invert1,invert2,geo1,geo2,collide,state1,state2);
        pair_energy_derivative(natom1,natom2,invert1,invert2,geo1,geo2,collide,state1,state2,grad1_minus,grad2_minus,dummy,dummy,dummy);
        state2[i] = state0;

        grad2[i] = (energy_plus - energy_minus)/(2.0*STEP);
        for(int j=0 ; j<7 ; j++)
        { hess22[i+j*7] = (grad2_plus[j] - grad2_minus[j])/(2.0*STEP); }
    }
}

// calculate a lattice-aligned bounding box for a molecule
void bound_box(int natom, // number of atoms in the molecule
               int invert, // inversion of the molecule {+1, -1}
               double *geo, // atomic coordinates of the molecule [3*natom]
               double *state, // state of the molecule [7]
               double *reclat, // reciprocal lattice vectors [6]
               double *box) // bounding box [6]
{
    box[0] = box[2] = box[4] = INFINITY;
    box[1] = box[3] = box[5] = -INFINITY;

    for(int i=0 ; i<natom ; i++)
    {
        double local[3], global[3], latpos[3];
        local[0] = invert*geo[3*i]; local[1] = invert*geo[1+3*i]; local[2] = invert*geo[2+3*i];
        position(local, state, global);
        latpos[0] = reclat[0]*global[0] + reclat[1]*global[1] + reclat[2]*global[2];
        latpos[1] = reclat[3]*global[1] + reclat[4]*global[2];
        latpos[2] = reclat[5]*global[2];

        if(latpos[0] < box[0]) { box[0] = latpos[0]; }
        if(latpos[0] > box[1]) { box[1] = latpos[0]; }
        if(latpos[1] < box[2]) { box[2] = latpos[1]; }
        if(latpos[1] > box[3]) { box[3] = latpos[1]; }
        if(latpos[2] < box[4]) { box[4] = latpos[2]; }
        if(latpos[2] > box[5]) { box[5] = latpos[2]; }
    }
}

// form reciprocal lattice vectors
void reciprocal(double *latvec, // lattice vectors (1st: x, 2nd: x-y, 3rd: x-y-z) [6]
                double *reclat) // reciprocal lattice vectors (1st: x-y-z, 2nd: y-z, 3rd: z) [6]
{
    double wt = 1.0/(latvec[0]*latvec[2]*latvec[5]);
    reclat[0] = wt*latvec[2]*latvec[5];
    reclat[1] = -wt*latvec[1]*latvec[5];
    reclat[2] = wt*(latvec[1]*latvec[4] - latvec[2]*latvec[3]);
    reclat[3] = wt*latvec[5]*latvec[0];
    reclat[4] = -wt*latvec[4]*latvec[0];
    reclat[5] = wt*latvec[0]*latvec[2];
}

// energy function that we are minimizing to relax the molecular crystal
double total_energy(struct molecular_crystal *xtl, // description of the crystal being optimized
                    double *state) // the crystal's state vector [6+7*xtl->nmol]
{
    double energy = fabs(state[0]*state[2]*state[5]);

    // construct reciprocal lattice vectors (1st: x-y-z, 2nd: y-z, 3rd: z)
    double reclat[6];
    reciprocal(state, reclat);

    // calculate buffer (lattice-aligned bounding box for the interaction sphere)
    double buffer[3];
    buffer[0] = INTERACTION_CUTOFF/sqrt(state[0]*state[0]);
    buffer[1] = INTERACTION_CUTOFF/sqrt(state[1]*state[1] + state[2]*state[2]);
    buffer[2] = INTERACTION_CUTOFF/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);

    // loop over molecules in the central unit cell
    for(int i=0 ; i<xtl->nmol ; i++)
    {
        double *state1 = state + 6 + 7*i;

        // calculate lattice-aligned bounding box for molecule 1
        double box1[6];
        bound_box(xtl->natom[xtl->type[i]], xtl->invert[i], xtl->geometry[xtl->type[i]], state1, reclat, box1);

        // loop over molecules in 2nd unit cell
        for(int j=i ; j<xtl->nmol ; j++)
        {
            // calculate lattice-aligned bounding box for molecule 2
            double box2[6];
            bound_box(xtl->natom[xtl->type[j]], xtl->invert[j], xtl->geometry[xtl->type[j]], state + 6 + 7*j, reclat, box2);

            // adjust lattice vector summation for pair of molecules
            int latmin1, latmin2, latmin3, latmax1, latmax2, latmax3;
            latmin1 = floor(box1[0] - box2[1] - buffer[0]);
            latmin2 = floor(box1[2] - box2[3] - buffer[1]);
            latmin3 = floor(box1[4] - box2[5] - buffer[2]);
            latmax1 = ceil(box1[1] - box2[0] + buffer[0]);
            latmax2 = ceil(box1[3] - box2[2] + buffer[1]);
            latmax3 = ceil(box1[5] - box2[4] + buffer[2]);

	    //printf("%d %d %d %d %d %d", latmin1, latmin2,latmin3,latmax1,latmax2,latmax3);
            // test if a crystal is too packed to continue
            if( (latmax1-latmin1+1)*(latmax1-latmin1+1)*(latmax1-latmin1+1) > MAX_CELL_SUM)
            { return INFINITY; }

            // loop over interacting unit cells
            for(int k=latmin1 ; k<=latmax1 ; k++)
            for(int l=latmin2 ; l<=latmax2 ; l++)
            for(int m=latmin3 ; m<=latmax3 ; m++)
            {
                // molecules don't interact with themselves & only consider unique pairs in central cell
                if(i == j && k == 0 && l == 0 && m == 0)
                { continue; }

                // reduced weights in central cell & between periodic images
                double wt = 2.0;
                if(i == j || (k == 0 && l == 0 && m == 0))
                { wt = 1.0; }

                // shift center of the 2nd molecule
                double state2[7];
                for(int n=0 ; n<7 ; n++)
                { state2[n] = state[n+6+7*j]; }
                state2[0] += k*state[0] + l*state[1] + m*state[3];
                state2[1] += l*state[2] + m*state[4];
                state2[2] += m*state[5];

		// Compute the distance between mol pairs
		double pair_dist = (state2[0] - state1[0]) * (state2[0] - state1[0]) 
		    + (state2[1] - state1[1]) * (state2[1] - state1[1])
		    + (state2[2] - state1[2]) * (state2[2] - state1[2]);
		pair_dist = sqrt(pair_dist);
		    
		// Pairs are too far apart - interaction is impossible 
		if(pair_dist > INTERACTION_CUTOFF + xtl->mol_length[i] + xtl->mol_length[j])
		{ continue; }

                energy += wt*pair_energy(xtl->natom[xtl->type[i]], xtl->natom[xtl->type[j]],
                                      xtl->invert[i], xtl->invert[j],
                                      xtl->geometry[xtl->type[i]], xtl->geometry[xtl->type[j]],
                                      xtl->collide[xtl->type[i]][xtl->type[j]],
                                      state1, state2);
                if(energy == INFINITY)
                { return energy; }
            }
        }
    }
    return energy;
}
void total_energy_derivative(struct molecular_crystal *xtl, // description of the crystal being optimized
                             double *state, // the crystal's state vector [6+7*xtl->nmol]
                             double *grad, // 1st derivative of the total energy [6+7*xtl->nmol]
                             double *hess) // 2nd derivative of the total energy [(6+7*xtl->nmol)*(6+7*xtl->nmol)]
{
    int size = 6+7*xtl->nmol;
    double sign0 = 1.0, sign2 = 1.0, sign5 = 1.0;
    if(state[0] < 0.0)
    { sign0 = -1.0; }
    if(state[2] < 0.0)
    { sign2 = -1.0; }
    if(state[5] < 0.0)
    { sign5 = -1.0; }

    for(int i=0 ; i<size ; i++)
    { grad[i] = 0.0; }
    grad[0] = fabs(state[2]*state[5])*sign0;
    grad[2] = fabs(state[0]*state[5])*sign2;
    grad[5] = fabs(state[0]*state[2])*sign5;

    for(int i=0 ; i<size*size ; i++)
    { hess[i] = 0.0; }
    hess[0 + 2*size] = fabs(state[5])*sign0*sign2;
    hess[2 + 0*size] = fabs(state[5])*sign0*sign2;
    hess[0 + 5*size] = fabs(state[2])*sign0*sign5;
    hess[5 + 0*size] = fabs(state[2])*sign0*sign5;
    hess[2 + 5*size] = fabs(state[0])*sign2*sign5;
    hess[5 + 2*size] = fabs(state[0])*sign2*sign5;

    // construct reciprocal lattice vectors (1st: x-y-z, 2nd: y-z, 3rd: z)
    double reclat[6];
    reciprocal(state, reclat);

    // calculate buffer (lattice-aligned bounding box for the interaction sphere)
    double buffer[3];
    buffer[0] = INTERACTION_CUTOFF/sqrt(state[0]*state[0]);
    buffer[1] = INTERACTION_CUTOFF/sqrt(state[1]*state[1] + state[2]*state[2]);
    buffer[2] = INTERACTION_CUTOFF/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);

    // loop over molecules in the central unit cell
    for(int i=0 ; i<xtl->nmol ; i++)
    {
        double *state1 = state + 6 + 7*i;

        // calculate lattice-aligned bounding box for molecule 1
        double box1[6];
        bound_box(xtl->natom[xtl->type[i]], xtl->invert[i], xtl->geometry[xtl->type[i]], state1, reclat, box1);

        // loop over molecules in 2nd unit cell
        for(int j=i ; j<xtl->nmol ; j++)
        {
            // calculate lattice-aligned bounding box for molecule 2
            double box2[6];
            bound_box(xtl->natom[xtl->type[j]], xtl->invert[j], xtl->geometry[xtl->type[j]], state + 6 + 7*j, reclat, box2);

            // adjust lattice vector summation for pair of molecules
            int latmin1, latmin2, latmin3, latmax1, latmax2, latmax3;
            latmin1 = floor(box1[0] - box2[1] - buffer[0]);
            latmin2 = floor(box1[2] - box2[3] - buffer[1]);
            latmin3 = floor(box1[4] - box2[5] - buffer[2]);
            latmax1 = ceil(box1[1] - box2[0] + buffer[0]);
            latmax2 = ceil(box1[3] - box2[2] + buffer[1]);
            latmax3 = ceil(box1[5] - box2[4] + buffer[2]);

            // test if a crystal is too packed to continue
            if( (latmax1-latmin1+1)*(latmax1-latmin1+1)*(latmax1-latmin1+1) > MAX_CELL_SUM)
            {
                for(int n=0 ; n<size ; n++)
                { grad[n] = NAN; }
                return;
            }

            // loop over interacting unit cells
            for(int k=latmin1 ; k<=latmax1 ; k++)
            for(int l=latmin2 ; l<=latmax2 ; l++)
            for(int m=latmin3 ; m<=latmax3 ; m++)
            {
                // molecules don't interact with themselves & only consider unique pairs in central cell
                if(i >= j && k == 0 && l == 0 && m == 0)
                { continue; }

                // reduced weights in central cell & between periodic images
                double wt = 2.0;
                if(i == j || (k == 0 && l == 0 && m == 0))
                { wt = 1.0; }

                // shift center of the 2nd molecule
                double state2[7];
                for(int n=0 ; n<7 ; n++)
                { state2[n] = state[n+6+7*j]; }
                state2[0] += k*state[0] + l*state[1] + m*state[3];
                state2[1] += l*state[2] + m*state[4];
                state2[2] += m*state[5];

		// Compute the distance between mol pairs
		double pair_dist = (state2[0] - state1[0]) * (state2[0] - state1[0]) +
		    (state2[1] - state1[1]) * (state2[1] - state1[1])+
		    (state2[2] - state1[2]) * (state2[2] - state1[2]);
		pair_dist = sqrt(pair_dist);
		// Pairs are too far apart - interaction is impossible 
		if(pair_dist > INTERACTION_CUTOFF + xtl->mol_length[i] + xtl->mol_length[j])
		{ continue; }

                double grad1[7], grad2[7], hess11[49], hess22[49], hess12[49];
                pair_energy_derivative(xtl->natom[xtl->type[i]], xtl->natom[xtl->type[j]],
                                       xtl->invert[i], xtl->invert[j],
                                       xtl->geometry[xtl->type[i]], xtl->geometry[xtl->type[j]],
                                       xtl->collide[xtl->type[i]][xtl->type[j]],
                                       state1, state2, grad1, grad2, hess11, hess22, hess12);
                for(int n=0 ; n<7 ; n++)
                { grad1[n] *= wt; grad2[n] *= wt; }
                for(int n=0 ; n<49 ; n++)
                { hess11[n] *= wt; hess12[n] *= wt; hess22[n] *= wt; }
                for(int n=0 ; n<7 ; n++)
                {
                    grad[n+6+i*7] += grad1[n];
                    grad[n+6+j*7] += grad2[n];
                }
                grad[0] += k*grad2[0];
                grad[1] += l*grad2[0];
                grad[2] += l*grad2[1];
                grad[3] += m*grad2[0];
                grad[4] += m*grad2[1];
                grad[5] += m*grad2[2];

                for(int n=0 ; n<7 ; n++)
                for(int o=0 ; o<7 ; o++)
                {
                    hess[n+6+i*7 + (o+6+i*7)*size] += hess11[n + o*7];
                    hess[n+6+j*7 + (o+6+j*7)*size] += hess22[n + o*7];
                    hess[n+6+i*7 + (o+6+j*7)*size] += hess12[n + o*7];
                    hess[n+6+j*7 + (o+6+i*7)*size] += hess12[o + n*7];
                }
                for(int n=0 ; n<7 ; n++)
                {
                    hess[0 + (n+6+i*7)*size] += k*hess12[n + 0*7];
                    hess[1 + (n+6+i*7)*size] += l*hess12[n + 0*7];
                    hess[2 + (n+6+i*7)*size] += l*hess12[n + 1*7];
                    hess[3 + (n+6+i*7)*size] += m*hess12[n + 0*7];
                    hess[4 + (n+6+i*7)*size] += m*hess12[n + 1*7];
                    hess[5 + (n+6+i*7)*size] += m*hess12[n + 2*7];
                    hess[0 + (n+6+j*7)*size] += k*hess22[n + 0*7];
                    hess[1 + (n+6+j*7)*size] += l*hess22[n + 0*7];
                    hess[2 + (n+6+j*7)*size] += l*hess22[n + 1*7];
                    hess[3 + (n+6+j*7)*size] += m*hess22[n + 0*7];
                    hess[4 + (n+6+j*7)*size] += m*hess22[n + 1*7];
                    hess[5 + (n+6+j*7)*size] += m*hess22[n + 2*7];

                    hess[n+6+i*7 + 0*size] += k*hess12[n + 0*7];
                    hess[n+6+i*7 + 1*size] += l*hess12[n + 0*7];
                    hess[n+6+i*7 + 2*size] += l*hess12[n + 1*7];
                    hess[n+6+i*7 + 3*size] += m*hess12[n + 0*7];
                    hess[n+6+i*7 + 4*size] += m*hess12[n + 1*7];
                    hess[n+6+i*7 + 5*size] += m*hess12[n + 2*7];
                    hess[n+6+j*7 + 0*size] += k*hess22[n + 0*7];
                    hess[n+6+j*7 + 1*size] += l*hess22[n + 0*7];
                    hess[n+6+j*7 + 2*size] += l*hess22[n + 1*7];
                    hess[n+6+j*7 + 3*size] += m*hess22[n + 0*7];
                    hess[n+6+j*7 + 4*size] += m*hess22[n + 1*7];
                    hess[n+6+j*7 + 5*size] += m*hess22[n + 2*7];
                }
                hess[0 + 0*size] += k*k*hess22[0 + 0*7];
                hess[1 + 0*size] += k*l*hess22[0 + 0*7];
                hess[2 + 0*size] += k*l*hess22[1 + 0*7];
                hess[3 + 0*size] += k*m*hess22[0 + 0*7];
                hess[4 + 0*size] += k*m*hess22[1 + 0*7];
                hess[5 + 0*size] += k*m*hess22[2 + 0*7];
                hess[0 + 1*size] += l*k*hess22[0 + 0*7];
                hess[1 + 1*size] += l*l*hess22[0 + 0*7];
                hess[2 + 1*size] += l*l*hess22[1 + 0*7];
                hess[3 + 1*size] += l*m*hess22[0 + 0*7];
                hess[4 + 1*size] += l*m*hess22[1 + 0*7];
                hess[5 + 1*size] += l*m*hess22[2 + 0*7];
                hess[0 + 2*size] += l*k*hess22[0 + 1*7];
                hess[1 + 2*size] += l*l*hess22[0 + 1*7];
                hess[2 + 2*size] += l*l*hess22[1 + 1*7];
                hess[3 + 2*size] += l*m*hess22[0 + 1*7];
                hess[4 + 2*size] += l*m*hess22[1 + 1*7];
                hess[5 + 2*size] += l*m*hess22[2 + 1*7];
                hess[0 + 3*size] += m*k*hess22[0 + 0*7];
                hess[1 + 3*size] += m*l*hess22[0 + 0*7];
                hess[2 + 3*size] += m*l*hess22[1 + 0*7];
                hess[3 + 3*size] += m*m*hess22[0 + 0*7];
                hess[4 + 3*size] += m*m*hess22[1 + 0*7];
                hess[5 + 3*size] += m*m*hess22[2 + 0*7];
                hess[0 + 4*size] += m*k*hess22[0 + 1*7];
                hess[1 + 4*size] += m*l*hess22[0 + 1*7];
                hess[2 + 4*size] += m*l*hess22[1 + 1*7];
                hess[3 + 4*size] += m*m*hess22[0 + 1*7];
                hess[4 + 4*size] += m*m*hess22[1 + 1*7];
                hess[5 + 4*size] += m*m*hess22[2 + 1*7];
                hess[0 + 5*size] += m*k*hess22[0 + 2*7];
                hess[1 + 5*size] += m*l*hess22[0 + 2*7];
                hess[2 + 5*size] += m*l*hess22[1 + 2*7];
                hess[3 + 5*size] += m*m*hess22[0 + 2*7];
                hess[4 + 5*size] += m*m*hess22[1 + 2*7];
                hess[5 + 5*size] += m*m*hess22[2 + 2*7];

                if(isnan(grad1[0]))
                {
                    for(int n=0 ; n<size ; n++)
                    { grad[n] = NAN; }
                    return;
                }
            }
        }
    }
}
// numerical derivatives for debugging purposes
void total_energy_derivative_test(struct molecular_crystal *xtl, // description of the crystal being optimized
                                  double *state, // the crystal's state vector [6+7*xtl->nmol]
                                  double *grad, // 1st derivative of the total energy [6+7*xtl->nmol]
                                  double *hess) // 2nd derivative of the total energy [(6+7*xtl->nmol)*(6+7*xtl->nmol)]
{
    int num = 6+7*xtl->nmol;
    double *grad_plus = (double*)malloc(sizeof(double)*num);
    double *grad_minus = (double*)malloc(sizeof(double)*num);
    double *dummy = (double*)malloc(sizeof(double)*num*num);
    for(int i=0 ; i<num ; i++)
    {
        double state0 = state[i];
        state[i] += STEP;
        double energy_plus = total_energy(xtl,state);
        total_energy_derivative(xtl,state,grad_plus,dummy);
        state[i] -= 2.0*STEP;
        double energy_minus = total_energy(xtl,state);
        total_energy_derivative(xtl,state,grad_minus,dummy);
        state[i] = state0;

        grad[i] = (energy_plus - energy_minus)/(2.0*STEP);

        for(int j=0 ; j<num ; j++)
        { hess[i+j*num] = (grad_plus[j] - grad_minus[j])/(2.0*STEP); }
    }
    free(grad_plus);
    free(grad_minus);
    free(dummy);
}

// renormalize quaternions, canonicalize lattice vectors, & recenter molecules
int renormalize(int nmol, // number of molecules in the state vector
                 double *state,
		 struct molecular_crystal *xtl) // state vector [6+7*nmol]
{
    // renormalize quaternions
    for(int i=0 ; i<nmol ; i++)
    {
        double renorm = 1.0/sqrt(state[9+7*i]*state[9+7*i] + state[10+7*i]*state[10+7*i]
                                + state[11+7*i]*state[11+7*i] + state[12+7*i]*state[12+7*i]);
        for(int j=0 ; j<4 ; j++)
        { state[9+j+7*i] *= renorm; }
    }

    // canonicalize lattice vectors
    int num;
    /*    if(state[2]*state[4] > 0.0)
    { num = state[4]/state[2] + 0.5; }
    else
    { num = state[4]/state[2] - 0.5; }
    state[3] -= num*state[1];
    state[4] -= num*state[2];

    if(state[0]*state[1] > 0.0)
    { num = state[1]/state[0] + 0.5; }
    else
    { num = state[1]/state[0] - 0.5; }
    state[1] -= num*state[0];

    if(state[0]*state[3] > 0.0)
    { num = state[3]/state[0] + 0.5; }
    else
    { num = state[3]/state[0] - 0.5; }
    state[3] -= num*state[0];
    
    // recenter molecules
    double reclat[6];
    reciprocal(state, reclat);
    for(int i=0 ; i<nmol ; i++)
    {
        double wt = state[6+7*i]*reclat[0] + state[7+7*i]*reclat[1] + state[8+7*i]*reclat[2];
        if(wt > 0.0) { num = wt + 0.5; }
        else { num = wt - 0.5; }
        state[6+7*i] -= num*state[0];

        wt = state[7+7*i]*reclat[3] + state[8+7*i]*reclat[4];
        if(wt > 0.0) { num = wt + 0.5; }
        else { num = wt - 0.5; }
        state[6+7*i] -= num*state[1];
        state[7+7*i] -= num*state[2];

        wt = state[8+7*i]*reclat[5];
        if(wt > 0.0) { num = wt + 0.5; }
        else { num = wt - 0.5; }
        state[6+7*i] -= num*state[3];
        state[7+7*i] -= num*state[4];
        state[8+7*i] -= num*state[5];
    }
    */

    // APPLY SYMMETRY OPERATIONS TO THE STATE VECTOR HERE (SYMMETRY INFO MUST BE INJECTED TO THIS POINT)
    if(xtl->spg != 0)
    {
	//symmetrize_state(state, xtl->invert, xtl->nmol, xtl->spg);

#ifdef ROPT_DEBUG
	crystal xtal;
	int N = xtl->natom[0];
	int Z = xtl->nmol;
	int size = 6 + 7 * xtl->nmol;
	xtal.num_atoms_in_molecule = N;
	xtal.Z = Z;
	xtal.Xcord = malloc(sizeof(float) * Z * N);
	xtal.Ycord = malloc(sizeof(float) * Z * N);
	xtal.Zcord = malloc(sizeof(float) * Z * N);
	xtal.atoms = malloc(sizeof(char) * Z * N * 2);
	for(int i = 0; i < Z * N; i++)
	{
	    xtal.atoms[2*i    ] = 'C';
	    xtal.atoms[2*i + 1] = ' ';
	}

	// Before symmetrization
	state_2_xtal(&xtal, state, xtl);
	//bring_all_molecules_to_first_cell(&xtal);
	//print_crystal(&xtal);
	xtal.num_atoms_in_molecule = N;
	//print_crystal(&xtal);
	int spg = detect_spg_using_spglib(&xtal);
	double temp[size];
	for(int i = 0; i < size; i++)
	    temp[i] = state[i];

	symmetrize_state(state, xtl->invert, xtl->nmol, xtl->spg);

	// After symmetrization
	state_2_xtal(&xtal, state, xtl);
	int spg2 = detect_spg_using_spglib(&xtal);
	printf("spg before = %d, spg after = %d\n", spg, spg2);

	// Check if symmetrize_state() is stationary
	/*	
	symmetrize_state(state, xtl->invert, xtl->nmol, xtl->spg);
	for(int i = 0; i < size; i++)
	{
	    if(fabs(state[i] - temp[i]) > 0.1)
	    {
		printf("Error in symmetrize_state: not stationary\n");
		print_state(state, size);
		print_state(temp, size);
		exit(1);
	    }	    
	}
	*/

	if(spg2 < spg)
	{
	   print_crystal(&xtal);
	   printf("Symmetrization failed.");
	   exit(1);
	}

	if(xtl->spg != spg2)
	{ return 1; }

	free(xtal.Xcord);
	free(xtal.Ycord);
	free(xtal.Zcord);
	free(xtal.atoms);
#endif
    }
    return 0;
}

// objective function to optimize the volume of the molecular crystal
double volume_search(double x, // optimization variable [0,1]
                     struct molecular_crystal *xtl, // description of the crystal being optimized
                     double *state, // the crystal's state vector [6+7*xtl->nmol]
                     double *min, // smallest scale factor to be considered [1]
                     double *max, // largest scale factor to be considered [1]
                     double *dummy, // dummy variable to match argument list w/ other objective function
                     double *dummy2, // dummy variable to match argument list w/ other objective function
                     double *work) // work vector [6+7*xtl->nmol]
{
    int size = 6+7*xtl->nmol;
    double scale0 = (1.0-x)*min[0] + x*max[0];
    for(int i=0 ; i<size ; i++)
    { work[i] = scale0*state[i]; }
    return total_energy(xtl, work);
}

// objective function for the Tikhonov-regularized line search
double quad_search(double x, // optimization variable [0,1]
                   struct molecular_crystal *xtl, // description of the crystal being optimized
                   double *state, // the crystal's state vector in normal coordinates [6+7*xtl->nmol]
                   double *tik_min, // minimum Tikhonov regularization parameter [1]
                   double *grad, // gradient in normal coordinates [6+7*xtl->nmol]
                   double *eval, // eigenvalues of the Hessian matrix [6+7*xtl->nmol]
                   double *evec, // eigenvectors of the Hessian matrix [(6+7*xtl->nmol)*(6+7*xtl->nmol)]
                   double *work) // work vector [13+14*xtl->nmol]
{
    int size = 6+7*xtl->nmol, sym_test;
    do
    {
        for(int i=0 ; i<size ; i++)
        {
            work[i] = state[i];
            work[i+size] = -x*grad[i]/(x*eval[i] + *tik_min);
        }
        char notrans = 'N';
        int inc = 1;
        double one = 1.0;
        dgemv_(&notrans, &size, &size, &one, evec, &size, work+size, &inc, &one, work, &inc);
        sym_test = renormalize(xtl->nmol, work, xtl);
        if(sym_test)
        {
            printf("symmetrization failed, reducing step size (%e -> %e)\n",x,0.5*x);
            x *= 0.5;
        }
    } while(sym_test);
    return total_energy(xtl, work);
}

// line optimizer for both objective functions (Golden-section search)
// code adapted from a Python implementation at https://en.wikipedia.org/wiki/Golden-section_search
double line_optimize(struct molecular_crystal *xtl, // description of the crystal being optimized
                     int nstep, // number of optimization steps
                     double *state, // the crystal's state vector in normal coordinates [6+7*xtl->nmol]
                     double (*fptr)(double, struct molecular_crystal*, double*, double*, double*, double*, double*, double*),
                     double *vec1,
                     double *vec2,
                     double *vec3,
                     double *vec4,
                     double *vec5)
{
    double invphi = (sqrt(5.0) - 1.0)*0.5, invphi2 = (3.0 - sqrt(5.0))*0.5;
    double a = 0.0, c = invphi2, d = invphi, h = 1.0;
    double yc = fptr(c, xtl, state, vec1, vec2, vec3, vec4, vec5);
    double yd = fptr(d, xtl, state, vec1, vec2, vec3, vec4, vec5);

    for(int i=0 ; i<nstep ; i++)
    {
        if(yc < yd)
        {
            d = c;
            yd = yc;
            h *= invphi;
            c = a + invphi2*h;
            yc = fptr(c, xtl, state, vec1, vec2, vec3, vec4, vec5);
        }
        else
        {
            a = c;
            c = d;
            yc = yd;
            h *= invphi;
            d = a + invphi*h;
            yd = fptr(d, xtl, state, vec1, vec2, vec3, vec4, vec5);
        }
    }

    // crappy hack to load the minimizer into the workspace (extra computation of objective function)
    double min = c, ymin = yc;
    if(yd < yc) { min = d; ymin = yd; }
    fptr(min, xtl, state, vec1, vec2, vec3, vec4, vec5);

    // extract minimizer from the workspace (a hack) & renormalize quaternions
    int size = 6+7*xtl->nmol;
    for(int i=0 ; i<size ; i++)
    { state[i] = vec5[i]; }
    return ymin;
}

// main loop of crystal optimization
Opt_status optimize(struct molecular_crystal *xtl, // description of the crystal being optimized
              double *state, // the crystal's state vector to be updated [6+7*xtl->nmol]
              Opt_settings set)
{
    int family = set.cell_family;
    int size = 6+7*xtl->nmol;
    double *workspace = (double*)malloc(sizeof(double)*size*2);

#ifdef ROPT_DEBUG
    clock_t start = clock();

    printf("\nStarted optimization with settings:\n");
    printf("Space group : %d\n", set.spg);
    printf("Cell family: %d\n", set.cell_family);
    printf("Max iteration: %d\n", set.max_iteration);
#endif
    
    // construct a constraint matrix for high-symmetry lattice vectors
    double constraint_mat[36];
    for(int i=0 ; i<36 ; i++)
    { constraint_mat[i] = 0.0; }
    switch(family)
    {
        case 1:
        constraint_mat[0 + 0*6] = constraint_mat[2 + 2*6] = constraint_mat[3 + 3*6] = constraint_mat[5 + 5*6] = 1.0;
        break;

        case 2:
        constraint_mat[0 + 0*6] = constraint_mat[2 + 2*6] = constraint_mat[5 + 5*6] = 1.0;
        break;

        case 3:
        constraint_mat[5 + 5*6] = 1.0;
        constraint_mat[0 + 0*6] = constraint_mat[2 + 0*6] = constraint_mat[0 + 2*6] = constraint_mat[2 + 2*6] = 0.5;
        break;

        case 4:
        constraint_mat[5 + 5*6] = 1.0;
        constraint_mat[0 + 0*6] = 0.5;
        constraint_mat[1 + 1*6] = 0.125;
        constraint_mat[2 + 2*6] = 0.375;
        constraint_mat[0 + 1*6] = constraint_mat[1 + 0*6] = -0.25;
        constraint_mat[0 + 2*6] = constraint_mat[2 + 0*6] = sqrt(3.0)*0.25;
        constraint_mat[1 + 2*6] = constraint_mat[2 + 1*6] = -sqrt(3.0)*0.125;
        break;

        case 5:
        constraint_mat[0 + 0*6] = constraint_mat[2 + 0*6] = constraint_mat[5 + 0*6] = 1.0/3.0;
        constraint_mat[0 + 2*6] = constraint_mat[2 + 2*6] = constraint_mat[5 + 2*6] = 1.0/3.0;
        constraint_mat[0 + 5*6] = constraint_mat[2 + 5*6] = constraint_mat[5 + 5*6] = 1.0/3.0;
        break;
    }

    // find an overpacked volume
    double scale_min = 1.0, scale_max = 1.0;
    double energy = total_energy(xtl, state), energy_min = energy, energy_max = energy;
    while(energy_min != INFINITY)
    {
        scale_min *= 0.5;
        energy_min = volume_search(0.0, xtl, state, &scale_min, &scale_max, NULL, NULL, workspace);
    }

    // find an underpacked volume
    while(energy_max == INFINITY)
    {
        scale_max *= 2.0;
        energy_max = volume_search(1.0, xtl, state, &scale_min, &scale_max, NULL, NULL, workspace);
        if(energy_max == INFINITY)
        { scale_min = scale_max; }
    }
    do
    {
        energy = energy_max;
        scale_max *= 2.0;
        energy_max = volume_search(1.0, xtl, state, &scale_min, &scale_max, NULL, NULL, workspace);
    } while(energy > energy_max);

#ifdef ROPT_DEBUG
    printf("Started preliminary volume optimization\n");
#endif
    // preliminary volume optimization
    energy = line_optimize(xtl, GOLDEN_STEPS, state, volume_search, &scale_min, &scale_max, NULL, NULL, workspace);
#ifdef ROPT_DEBUG
    printf("Completed preliminary volume optimization\n");
#endif
 
    // main optimization loop
    int niter = 0; // iteration counter
    int lwork = -1, info, inc = 1, six = 6;
    char jobz = 'V', uplo = 'U', notrans = 'N', trans = 'T';
    double work0, one = 1.0, zero = 0.0;
    double *grad = (double*)malloc(sizeof(double)*size);
    double *hess = (double*)malloc(sizeof(double)*size*size);
    double *ev = (double*)malloc(sizeof(double)*size);
    dsyev_(&jobz, &uplo, &size, NULL, &size, NULL, &work0, &lwork, &info);
    lwork = (int)work0;
    if(lwork < 6*size)
    { lwork = 6*size; }
    double *work = (double*)malloc(sizeof(double)*lwork);
    double new_energy = energy, energy_reduction_estimate;
    do
    {
        // save previous energy
        energy = new_energy;

        // expand the total energy to 2nd order
        total_energy_derivative(xtl, state, grad, hess);

        // apply constraints to gradient & Hessian
        if(family != 0)
        {
            dgemv_(&notrans, &six, &six, &one, constraint_mat, &six, grad, &inc, &zero, work, &inc);
            for(int i=0 ; i<6 ; i++)
            { grad[i] = work[i]; }

            dgemm_(&notrans, &notrans, &six, &size, &six, &one, constraint_mat, &six, hess, &size, &zero, work, &six);
            for(int i=0 ; i<6 ; i++)
            for(int j=0 ; j<size ; j++)
            { hess[i+j*size] = work[i+j*6]; }

            dgemm_(&notrans, &notrans, &size, &six, &six, &one, hess, &size, constraint_mat, &six, &zero, work, &size);
            for(int i=0 ; i<size ; i++)
            for(int j=0 ; j<6 ; j++)
            { hess[i+j*size] = work[i+j*size]; }
        }

        // transform into normal modes of the quadratic approximant
        dsyev_(&jobz, &uplo, &size, hess, &size, ev, work, &lwork, &info);
        for(int i=0 ; i<size ; i++)
        { work[i] = grad[i]; }
        dgemv_(&trans, &size, &size, &one, hess, &size, work, &inc, &zero, grad, &inc);

        // start with a sensible minimum Tikhonov regularization value
        double min_tik = 0.0;
        for(int i=0 ; i<size ; i++)
        { min_tik += grad[i]*grad[i]; }
        min_tik = sqrt(min_tik)/(MAX_LINE_STEP*(double)xtl->nmol) - ev[0];
        if(min_tik < 0.0)
        { min_tik = 0.0; }

        // calculate an estimate of energy reduction for testing false convergence
        energy_reduction_estimate = 0.0;
        for(int i=0 ; i<size ; i++)
        { energy_reduction_estimate -= grad[i]*grad[i]/(ev[i] + min_tik); }

        // increase the Tikhonov value until the 1st searched point lowers the energy
        while(quad_search((3.0 - sqrt(5.0))*0.5, xtl, state, &min_tik, grad, ev, hess, workspace) > energy)
        { min_tik *= 2.0; }

        // perform a Tikhonov-regularized line search
        new_energy = line_optimize(xtl, GOLDEN_STEPS, state, quad_search, &min_tik, grad, ev, hess, workspace);

        // strictly enforce constraints at the end of every optimization step
        switch(family)
        {
            case 1:
            state[1] = state[4] = 0.0;
            break;

            case 2:
            state[1] = state[3] = state[4] = 0.0;
            break;

            case 3:
            state[1] = state[3] = state[4] = 0.0;
            state[2] = state[0];
            break;

            case 4:
            state[3] = state[4] = 0.0;
            state[1] = -0.5*state[0];
            state[2] = sqrt(3.0)*0.5*state[0];
            break;

            case 5:
            state[1] = state[3] = state[4] = 0.0;
            state[2] = state[0];
            state[5] = state[0];
            break;
        }

	niter++;

#ifdef ROPT_DEBUG
	printf("Completed iteration #%5d: New energy =  %lf, Delta E = %lf\n",
	       niter, new_energy, (energy - new_energy));
	fflush(stdout);
#endif

	//}while(1);
    } while((energy - new_energy) > OPTIMIZATION_TOLERANCE*fabs(new_energy) &&
     niter < set.max_iteration);

#ifdef ROPT_DEBUG
    clock_t diff = clock() - start;
    printf("optimization time = %e s\n", (double)diff/(double)CLOCKS_PER_SEC);
#endif
    free(workspace);
    free(grad);
    free(hess);
    free(ev);
    free(work);

    if(niter == set.max_iteration)
    {  return ITER_LIMIT; }

    if(energy_reduction_estimate < -ENERGY_REDUCTION_TOLERANCE*xtl->nmol)
    { return FALSE_CONVERGENCE; }

    return SUCCESS;
    
}

// adapted from pseudocode in "Converting a Rotation Matrix to a Quaternion" by Mike Day
void matrix2quaternion(double *rot, double *quat)
{
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

// NOTE: rows/columns of the cutoff_matrix are over all atoms in the unit cell
Opt_status optimize_crystal(crystal *xtl, float *cutoff_matrix, int placeholder, Opt_settings set)
{
    set.cell_family = get_cell_type_from_spg(set.spg);
    int family = set.cell_family;
    // test of lattice vector formatting
    if(family != 0 && (xtl->lattice_vectors[0][1] != 0.0 ||
		       xtl->lattice_vectors[0][2] != 0.0 ||
		       xtl->lattice_vectors[1][2] != 0.0))
    {return MISC_FAILURE;}

    // allocate memory for the temporary crystal data structure & state vector
    struct molecular_crystal xtl2;
    xtl2.spg = set.spg;
    xtl2.ntype = 1;
    xtl2.nmol = xtl->Z;
    xtl2.natom = (int*)malloc(sizeof(int)*1);
    xtl2.natom[0] = xtl->num_atoms_in_molecule;
    xtl2.geometry = (double**)malloc(sizeof(double*)*1);
    xtl2.geometry[0] = (double*)malloc(sizeof(double)*3*xtl2.natom[0]);
    xtl2.collide = (double***)malloc(sizeof(double**)*1);
    xtl2.collide[0] = (double**)malloc(sizeof(double*)*1);
    xtl2.collide[0][0] = (double*)malloc(sizeof(double)*xtl2.natom[0]*xtl2.natom[0]);
    for(int i=0 ; i<xtl2.natom[0] ; i++)
    for(int j=0 ; j<xtl2.natom[0] ; j++)
    { xtl2.collide[0][0][j+i*xtl2.natom[0]] = cutoff_matrix[j+i*xtl2.natom[0]*xtl2.nmol]; }
    xtl2.type = (int*)malloc(sizeof(int)*xtl2.nmol);
    xtl2.invert = (int*)malloc(sizeof(int)*xtl2.nmol);
    xtl2.mol_length = (double *)malloc(sizeof(double)*xtl2.nmol);
    for(int i=0 ; i<xtl2.nmol ; i++)
    { xtl2.type[i] = 0; xtl2.invert[i] = 1; }
    double *state = (double*)malloc(sizeof(double)*(6 + 7*xtl2.nmol));

    // Compute molecule lengths
    for(int i = 0; i < xtl->Z; i++)
    {
	double cog[3] = {0};
	for(int j = 0; j < xtl->num_atoms_in_molecule; j++)
	{
	    cog[0] += xtl->Xcord[j];
	    cog[1] += xtl->Ycord[j];
	    cog[2] += xtl->Zcord[j];
	}
	
	for(int j =0; j < 3; j++)
	{ cog[j] /= xtl->num_atoms_in_molecule; }

	double max_dist = 0;
	double dist;
	for(int j = 0; j < xtl->num_atoms_in_molecule; j++)
	{
	    dist = sqrt( (xtl->Xcord[j] - cog[0])*(xtl->Xcord[j] - cog[0])
			+(xtl->Ycord[j] - cog[1])*(xtl->Ycord[j] - cog[1])
			+(xtl->Zcord[j] - cog[2])*(xtl->Zcord[j] - cog[2])
			);
	    if(dist > max_dist)
	    { max_dist = dist;}
	}
	xtl2.mol_length[i] = max_dist;
	//printf("mol length = %lf", max_dist);
    }
    
    // form rotation matrix to align lattice vectors (QR decomposition)
    double latvec[9], latvec2[9];
    for(int i=0 ; i<3 ; i++)
    for(int j=0 ; j<3 ; j++)
    { latvec[j + i*3] = xtl->lattice_vectors[i][j]; }
    for(int i=0 ; i<9 ; i++)
    { latvec2[i] = latvec[i]; }
    double tau[3];
    int lwork = xtl2.natom[0]+100, info, three = 3, one = 1;
    double *work = (double*)malloc(sizeof(double)*lwork);
    if(family == 0)
    { dgeqrf_(&three, &three, latvec2, &three, tau, work, &lwork, &info); }
    state[0] = latvec2[0 + 0*3];
    state[1] = latvec2[0 + 1*3];
    state[2] = latvec2[1 + 1*3];
    state[3] = latvec2[0 + 2*3];
    state[4] = latvec2[1 + 2*3];
    state[5] = latvec2[2 + 2*3];

    // isolate reference molecule
    char trans = 'T', notrans = 'N', left = 'L';
    for(int i=0 ; i<xtl2.natom[0] ; i++)
    {
        xtl2.geometry[0][0+3*i] = xtl->Xcord[i];
        xtl2.geometry[0][1+3*i] = xtl->Ycord[i];
        xtl2.geometry[0][2+3*i] = xtl->Zcord[i];
    }
    if(family == 0)
    { dormqr_(&left, &trans, &three, xtl2.natom, &three, latvec2, &three, tau, xtl2.geometry[0], &three, work, &lwork, &info); }

    // extract geometric center & center reference geometry
    double coord[3];
    coord[0] = coord[1] = coord[2] = 0.0;
    for(int i=0 ; i<xtl2.natom[0] ; i++)
    for(int j=0 ; j<3 ; j++)
    { coord[j] += xtl2.geometry[0][j+3*i]; }
    coord[0] /= (double)xtl2.natom[0];
    coord[1] /= (double)xtl2.natom[0];
    coord[2] /= (double)xtl2.natom[0];
    state[6] = coord[0];
    state[7] = coord[1];
    state[8] = coord[2];
    state[9] = 1.0;
    state[10] = 0.0;
    state[11] = 0.0;
    state[12] = 0.0;
    for(int i=0 ; i<xtl2.natom[0] ; i++)
    for(int j=0 ; j<3 ; j++)
    { xtl2.geometry[0][j+3*i] -= coord[j]; }

    // align & center the remaining molecules
    double *geo = (double*)malloc(sizeof(double)*xtl2.natom[0]*3);
    for(int i=1 ; i<xtl2.nmol ; i++)
    {
        // transform molecular coordinates to be consistent w/ lattice vectors
        for(int j=0 ; j<xtl2.natom[0] ; j++)
        {
            geo[0+j*3] = xtl->Xcord[j+i*xtl2.natom[0]];
            geo[1+j*3] = xtl->Ycord[j+i*xtl2.natom[0]];
            geo[2+j*3] = xtl->Zcord[j+i*xtl2.natom[0]];
        }
        if(family == 0)
        { dormqr_(&left, &trans, &three, xtl2.natom, &three, latvec2, &three, tau, geo, &three, work, &lwork, &info); }

        // extract center of molecule
        for(int j=0 ; j<3 ; j++)
        {
            state[j+6 + 7*i] = 0.0;
            for(int k=0 ; k<xtl2.natom[0] ; k++)
            { state[j+6 + 7*i] += geo[j+k*3]; }
            state[j+6 + 7*i] /= (double)xtl2.natom[0];
        }
        for(int j=0 ; j<xtl2.natom[0] ; j++)
        for(int k=0 ; k<3 ; k++)
        { geo[k + j*3] -= state[k+6 + 7*i]; }

        // form rotation matrix (Kabsch algorithm)
        char all = 'A';
        double sv[3], rot[9], U[9], VT[9], zero = 0.0, real_one = 1.0;
        dgemm_(&notrans, &trans, &three, &three, xtl2.natom, &real_one, geo, &three, xtl2.geometry[0], &three, &zero, rot, &three);
        dgesvd_(&all, &all, &three, &three, rot, &three, sv, U, &three, VT, &three, work, &lwork, &info);
        dgemm_(&notrans, &notrans, &three, &three, &three, &real_one, U, &three, VT, &three, &zero, rot, &three);

        // check for inversion (determinant test)
        double det = rot[0+0*3]*rot[1+1*3]*rot[2+2*3]
                    -rot[0+0*3]*rot[1+2*3]*rot[2+1*3]
                    +rot[0+1*3]*rot[1+2*3]*rot[2+0*3]
                    -rot[0+1*3]*rot[1+0*3]*rot[2+2*3]
                    +rot[0+2*3]*rot[1+0*3]*rot[2+1*3]
                    -rot[0+2*3]*rot[1+1*3]*rot[2+0*3];
        if(det < 0.0)
        {
            xtl2.invert[i] = -1;
            for(int j=0 ; j<9 ; j++)
            { rot[j] *= (double)xtl2.invert[i]; }
        }

        // extract quaternion vector
        matrix2quaternion(rot, state + 9 + 7*i);
    }

   // run the optimizer
    Opt_status status = optimize(&xtl2, state, set);

    // convert the result back to the original format
    for(int i=0 ; i<9 ; i++)
    { latvec[i] = 0.0; }
    latvec[0 + 0*3] = state[0];
    latvec[0 + 1*3] = state[1];
    latvec[1 + 1*3] = state[2];
    latvec[0 + 2*3] = state[3];
    latvec[1 + 2*3] = state[4];
    latvec[2 + 2*3] = state[5];
    if(family == 0)
    { dormqr_(&left, &notrans, &three, &three, &three, latvec2, &three, tau, latvec, &three, work, &lwork, &info); }

    for(int i=0 ; i<3 ; i++)
    for(int j=0 ; j<3 ; j++)
    { xtl->lattice_vectors[i][j] = latvec[j + i*3]; }

    for(int i=0 ; i<xtl2.nmol ; i++)
    for(int j=0 ; j<xtl2.natom[0] ; j++)
    {
        double local[3];
        for(int k=0 ; k<3 ; k++)
        { local[k] = (double)xtl2.invert[i]*xtl2.geometry[0][k+3*j]; }
        position(local, state + 6 + 7*i, coord);
        if(family == 0)
        { dormqr_(&left, &notrans, &three, &one, &three, latvec2, &three, tau, coord, &three, work, &lwork, &info); }

        xtl->Xcord[j+i*xtl2.natom[0]] = coord[0];
        xtl->Ycord[j+i*xtl2.natom[0]] = coord[1];
        xtl->Zcord[j+i*xtl2.natom[0]] = coord[2];
    }

    free_molecular_crystal(&xtl2);
    free(state);
    free(geo);
    free(work);

    return status;
}

void state_2_xtal(crystal *xtl, double *state, struct molecular_crystal *xtl2)
{
    double tau[3];
    char notrans = 'N', left = 'L';
    int lwork = 1000, info, three = 3, one = 1;
    double work;
    
    double latvec[9] = {0};
    latvec[0 + 0*3] = state[0];
    latvec[0 + 1*3] = state[1];
    latvec[1 + 1*3] = state[2];
    latvec[0 + 2*3] = state[3];
    latvec[1 + 2*3] = state[4];
    latvec[2 + 2*3] = state[5];
   
    //if(family == 0)
    //{ dormqr_(&left, &notrans, &three, &three, &three, latvec2, &three, tau, latvec, &three, &work, &lwork, &info); }

    for(int i=0 ; i<3 ; i++)
    for(int j=0 ; j<3 ; j++)
    { xtl->lattice_vectors[i][j] = latvec[j + i*3]; }

    int jmol = 0, imol = 0;
    double coord[3];
    for(int i=0 ; i<xtl2->nmol ; i++)
    {
        imol = xtl2->type[i];
        for(int j=0 ; j<xtl2->natom[imol] ; j++)
        {
            double local[3];
            for(int k=0 ; k<3 ; k++)
            { local[k] = (double)xtl2->invert[i]*xtl2->geometry[imol][k+3*j]; }
            position(local, state + 6 + 7*i, coord);
            //if(family == 0)
            //{ dormqr_(&left, &notrans, &three, &one, &three, latvec2, &three, tau, coord, &three, work, &lwork, &info); }

            xtl->Xcord[j+jmol] = coord[0];
            xtl->Ycord[j+jmol] = coord[1];
            xtl->Zcord[j+jmol] = coord[2];
        }
        jmol += xtl2->natom[imol];
    }
}

static int get_cell_type_from_spg(int spg)
{
    if (spg <= 2)
        return 0;
    else if (spg <= 15)
        return 1;
    else if (spg <= 74)
        return 2;
    else if (spg <= 142)
        return 3;
    else if (spg <= 194)
        return 4;
    else if (spg <= 230)
        return 5;

    return 0;
}
