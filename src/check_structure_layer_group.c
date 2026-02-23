#include <stdio.h>
#include <math.h>
#include <string.h>
#include "check_structure_layer_group.h"
#include "read_input.h"
#include "spg_generation.h"
#include "crystal_utils.h"
#include "combinatorics.h"
#include "molecule_placement.h"
#include "molecule_utils.h"
#include "algebra.h"
#include "spglib.h"

void check_layer_group(crystal* xtal)
{

    molecule *mol = (molecule*)malloc(sizeof(molecule));//store molecule
    read_geometry(mol);             //read molecule from geometry.in

    //recenter molecule to origin
    recenter_molecule(mol);
    int attempted_lg = xtal->spg;
    printf("The attempted layer group is %d\n",attempted_lg);
    int spglib_spg = detect_spg_using_spglib(xtal);
    printf( "#SPGLIB_detected_spacegroup = %d\n", spglib_spg);
    int spglib_lg[5];
    int match = 0;
    for (int i =0;i<5;i++)
    {
        spglib_lg[i] = spg_to_lg[spglib_spg-1][i];
        if (spglib_lg[i] == attempted_lg)
        {
            match == 1;
	    printf("match = 1\n");
        }
        printf("the lg is %d\n",spglib_lg[i]);
    }

   

    //if result in a higher symmetry
    if (match == 0)
    {
	for (int j = 0; j<5;j++)
	    {
		lg = spglib_lg[j];
		if (lg != 0)
		    {
		        crystal *tmp_xtal = (crystal *)malloc(sizeof(crystal) );
    			allocate_xtal(tmp_xtal,Z,N);
		        create_xtal_in_designated_lg(lg,xtal,tmp_xtal,mol);	
		


			free_xtal(tmp_xtal);		
		    }

	    }

    }


    //free_xtal(tmp_xtal);
}

void create_xtal_in_designated_lg(int lg,crystal* xtal,crystal* tmp_xtali,molecule* mol)
{    //create a temp xtal
    float Z = xtal->Z;
    float N = xtal->num_atoms_in_molecule;
    //crystal *tmp_xtal = (crystal *)malloc(sizeof(crystal) );
    //allocate_xtal(tmp_xtal,Z,N);
    
    float lattice_vec_a[3];
    float lattice_vec_b[3];
    float lattice_vec_c[3];
    float norm_a = 0;
    float norm_b = 0;
    float norm_c = 0;
    float mol_Xfrac[N]; //stores fractional coordinates of first mol
    float mol_Yfrac[N];
    float mol_Zfrac[N];
    // get a,b,and c and assign lattice vectors to tmp_xtal
    for (int i = 0; i < 3; i++)
	{
            tmp_xtal->lattice_vectors[0][i] = xtal->lattice_vectors[0][i];
	    tmp_xtal->lattice_vectors[1][i] = xtal->lattice_vectors[1][i];
	    tmp_xtal->lattice_vectors[2][i] = xtal->lattice_vectors[2][i];
	    norm_a = norm_a + pow(tmp_xtal->lattice_vectors[0][i],2);
	    norm_b = norm_b + pow(tmp_xtal->lattice_vectors[1][i],2);
	    norm_c = norm_c + pow(tmp_xtal->lattice_vectors[2][i],2);

	}
	norm_a = sqrt(norm_a);
        norm_b = sqrt(norm_b);
	norm_c = sqrt(norm_c);
   


    // store position of first mol in array, with frac coords
    for(int i = 0;i < N;i++)
	{
	    mol_Xfrac[i] = xtal->Xcord[i]/norm_a;
	    mol_Yfrac[i] = xtal->Ycord[i]/norm_b;
	    mol_Zfrac[i] = xtal->Zcord[i]/norm_c;
	}
 
    int hall_number = lg;
    apply_all_lg_symmetry_ops(tmp_xtal,
                              mol,
                              mol_Xfrac,
                              mol_Yfrac,
                              mol_Zfrac,
                              N,
                              hall_number);



}

/*


void main()
{
    int spglib_lg[5];
    for (int i =0;i<5;i++)
    {
        spglib_lg[i] = spg_to_lg[192][i];
        printf("the lg is %d\n",spglib_lg[i]);
    }
    

}
*/
