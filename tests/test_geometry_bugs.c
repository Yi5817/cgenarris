#include <math.h>
#include <stdio.h>
#include <string.h>

#include "algebra.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "lattice_generator.h"
#include "molecule_utils.h"
#include "spg_generation.h"
#include "combinatorics.h"

extern float TOL;
static int failures;

static void expect(int condition, const char *name)
{
    if(!condition)
    {
        fprintf(stderr, "FAIL: %s\n", name);
        failures++;
    }
}

/* Independent exhaustive image enumeration for these bounded fixtures. */
static int reference_accepts(const crystal *c, const float *cut)
{
    int n = c->Z * c->num_atoms_in_molecule;
    for(int a = 0; a < n; a++)
    for(int b = 0; b < n; b++)
    for(int i = -12; i <= 12; i++)
    for(int j = -12; j <= 12; j++)
    for(int k = -2; k <= 2; k++)
    {
        if(i == 0 && j == 0 && k == 0 &&
           a / c->num_atoms_in_molecule == b / c->num_atoms_in_molecule)
            continue;
        double dx = c->Xcord[a] - c->Xcord[b] -
            i*c->lattice_vectors[0][0] - j*c->lattice_vectors[1][0] -
            k*c->lattice_vectors[2][0];
        double dy = c->Ycord[a] - c->Ycord[b] -
            j*c->lattice_vectors[1][1] - k*c->lattice_vectors[2][1];
        double dz = c->Zcord[a] - c->Zcord[b] - k*c->lattice_vectors[2][2];
        if(dx*dx + dy*dy + dz*dz < (double)cut[a*n+b]*cut[a*n+b])
            return 0;
    }
    return 1;
}

static void test_contacts(void)
{
    float x[4] = {-.37f, .37f, -.37f, .37f};
    float y[4] = {0, 0, 8, 8}, z[4] = {0};
    char atoms[8] = {'H',' ','H',' ','H',' ','H',' '};
    float cut[16];
    for(int i = 0; i < 16; i++) cut[i] = 2.2f;
    crystal c = {.lattice_vectors={{10,0,0},{0,10,0},{0,0,10}},
        .Xcord=x, .Ycord=y, .Zcord=z, .atoms=atoms,
        .Z=2, .num_atoms_in_molecule=2};
    expect(!reference_accepts(&c, cut), "reference H2 clash");
    expect(!check_structure_with_vdw_matrix(c, cut, 4, 4), "H2 image clash");
    y[2] = y[3] = 5;
    expect(reference_accepts(&c, cut), "reference separated H2");
    expect(check_structure_with_vdw_matrix(c, cut, 4, 4), "separated H2");

    x[0]=-.7f; x[1]=.7f; x[2]=3.05f; x[3]=4.45f;
    memset(y, 0, sizeof(y)); c.lattice_vectors[0][0]=7.5f;
    for(int i = 0; i < 16; i++) cut[i] = 3.4f;
    expect(!reference_accepts(&c, cut), "reference carbon clash");
    expect(!check_structure_with_vdw_matrix(c, cut, 4, 4), "carbon cutoff");

    c.lattice_vectors[0][0]=10; c.lattice_vectors[1][0]=5.1961524f;
    c.lattice_vectors[1][1]=3; c.lattice_vectors[2][0]=115.47005f;
    c.lattice_vectors[2][1]=200; c.lattice_vectors[2][2]=10;
    x[0]=-1; x[1]=1; x[2]=15.373067f; x[3]=17.373067f;
    y[2]=y[3]=21; z[2]=z[3]=1;
    for(int i = 0; i < 16; i++) cut[i] = 2.2f;
    expect(!reference_accepts(&c, cut), "reference skew clash");
    expect(!check_structure_with_vdw_matrix(c, cut, 4, 4), "skew image clash");

    /* Zero molecular diameter must still check periodic self-images. */
    c.Z=1; c.num_atoms_in_molecule=1;
    float cell[3][3]={{3,0,0},{0,10,0},{0,0,10}};
    memcpy(c.lattice_vectors, cell, sizeof(cell)); cut[0]=3.4f;
    expect(!reference_accepts(&c, cut), "reference monatomic self clash");
    expect(!check_structure_with_vdw_matrix(c, cut, 1, 1), "self image cutoff");
}

static void test_species(void)
{
    float x[4]={-1,1,1,-1}, y[4]={0}, z[4]={0};
    char atoms[8]={'C',' ','O',' ','C',' ','O',' '};
    crystal c={.lattice_vectors={{10,0,0},{0,10,0},{0,0,10}},
        .Xcord=x,.Ycord=y,.Zcord=z,.atoms=atoms,.Z=2,.num_atoms_in_molecule=2};
    int overlap[2]={0,2};
    TOL=.1f;
    expect(!check_overlap_xtal(&c,overlap,2,2), "inverted CO periodic");
    expect(!check_overlap_xtal_cartesian(&c,overlap,2,2), "inverted CO Cartesian");
    atoms[4]='O'; atoms[6]='C';
    expect(check_overlap_xtal(&c,overlap,2,2), "valid species permutation");
    expect(check_overlap_xtal_cartesian(&c,overlap,2,2), "valid Cartesian permutation");
    combine_close_molecules(&c);
    expect(c.Z == 1 && fabsf(x[0]+1)<1e-6f && fabsf(x[1]-1)<1e-6f,
           "averaging preserves species correspondence");
    c.Z=2; atoms[0]='C'; atoms[2]='C'; atoms[3]='l';
    atoms[4]='C'; atoms[5]=' '; atoms[6]='C'; atoms[7]='l';
    x[0]=-1; x[1]=1; x[2]=1; x[3]=-1;
    expect(!check_overlap_xtal(&c,overlap,2,2), "both symbol characters matter");

    /* Each second atom is near the same first atom: no bijection exists. */
    for(int i=0; i<4; i++) { atoms[2*i]='C'; atoms[2*i+1]=' '; }
    x[0]=0; x[1]=1; x[2]=0; x[3]=.01f;
    expect(!check_overlap_xtal(&c,overlap,2,2), "reject many-to-one matching");
    /* A greedy assignment fails; a valid reassignment exists. */
    x[0]=0; x[1]=.15f; x[2]=.075f; x[3]=-.075f;
    expect(check_overlap_xtal(&c,overlap,2,2), "augmenting atom assignment");
}

static void test_alignment(void)
{
    float targets[][3]={{-1,0,0},{1,0,0},{1,.01f,0},{-1,.01f,0},{0,1,0}};
    for(int i=0; i<5; i++)
    {
        float a[3]={1,0,0}, b[3], rot[3][3], out[3];
        memcpy(b,targets[i],sizeof(b)); normalise_vector3(b);
        rotation_matrix_from_vectors(rot,a,b);
        vector3_mat3b3_multiply(rot,a,out);
        expect(cart_dist(out,b)<1e-5f, "vector alignment");
        expect(fabsf(det_mat3b3(rot)-1)<1e-5f, "proper rigid rotation");
    }
}

static void test_small_volumes(void)
{
    float lat[3][3];
    gen_orthorhombic_lattice(lat,8,.3f);
    expect(!isfinite(lat[0][0]), "impossible volume rejected");
    gen_orthorhombic_lattice(lat,27,.3f);
    expect(lat[0][0]==3 && lat[1][1]==3 && lat[2][2]==3,
           "boundary volume gives unique feasible diagonal");
    gen_orthorhombic_lattice(lat,27.0001f,.3f);
    expect(!isfinite(lat[0][0]) ||
           (lat[0][0]>=3 && lat[1][1]>=3 && lat[2][2]>=3),
           "near-boundary draw terminates safely");
}

int main(int argc, char **argv)
{
    (void)argv;
    if(argc > 1) test_small_volumes();
    else { test_contacts(); test_species(); test_alignment(); }
    return failures ? 1 : 0;
}
