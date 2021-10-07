#include <stdio.h>
#include <math.h>
#include <string.h>
#include "spglib.h"
#include "algebra.h"
#include "read_input.h"
#include "spg_generation.h"
#include "molecule_placement.h"
#include "molecule_utils.h"
#include "crystal_utils.h"
#include "randomgen.h"
#include "combinatorics.h"
#include "layer_group_position_database.h"

#define PI 3.141592653
extern unsigned int *seed2;


int cxtal_place_molecules(cocrystal *cxtal, Settings set, molecule *mol)
{
    float random_rot[3][3];
    float random_trans[3];
    int asym_atoms = cxtal->n_atoms/ set.Z;
    float frac[asym_atoms][3];
    float inv_lat_vec[3][3];
    inverse_mat3b3(inv_lat_vec, cxtal->lattice_vectors);
    mat3b3_transpose(inv_lat_vec, inv_lat_vec);

    // Get fractional coordinates for asym unit
    int at = 0;
    for(int m = 0; m < cxtal->n_mol_types; m++)
    {
        for(int st = 0; st < cxtal->stoic[m]; st++)
        {
            generate_random_rotation_matrix(random_rot);
            generate_random_translation_vector(random_trans);

            for(int mat = 0; mat < cxtal->n_atoms_in_mol[m]; mat++)
            {
                float atom[3] = {mol[m].X[mat], mol[m].Y[mat], mol[m].Z[mat]};
                vector3_mat3b3_multiply(random_rot, atom, atom);

                // Convert to fractional and translate
                vector3_mat3b3_multiply(inv_lat_vec, atom, atom);
                vector3_add(atom, random_trans, atom);

                memcpy(frac[at], atom, 3*sizeof(float));
                cxtal->atoms[2*at] = mol[m].atoms[2*mat];
                cxtal->atoms[2*at + 1] = mol[m].atoms[2*mat + 1];
                at++;
            }
        }
    }

    cxtal_apply_symmetry_ops(cxtal, frac);
    return 1;
}

/*
Apply symmetry operation to generate the full crystal.
The fractional coordinates of the asym unit is frac.
*/
void cxtal_apply_symmetry_ops(cocrystal *cxtal, float (*frac)[3])
{
    //get symmetry operations from spglib database
    double translations[192][3];
    int rotations[192][3][3];
    int hall_number = hall_number_from_spg(cxtal->spg);
    int num_of_operations = spg_get_symmetry_from_database(rotations,
                            translations, hall_number);
    int n_atoms_asym = cxtal->n_atoms / cxtal->Z;
    float lat_vec_trans[3][3];
    mat3b3_transpose(lat_vec_trans, cxtal->lattice_vectors);

    int at = 0;
    for(int op = 0; op < num_of_operations; op++)
    {
        float trans[3] = {translations[op][0],
			  translations[op][1],
			  translations[op][2]};
        int (*rot)[3] = rotations[op];

        for(int asym = 0; asym < n_atoms_asym; asym++)
        {
            float atom[3];
            memcpy(atom, frac[asym], 3*sizeof(float));
            vector3_intmat3b3_multiply(rot, atom, atom);
            vector3_add(trans, atom, atom);
            // Convert to cartesian
            vector3_mat3b3_multiply(lat_vec_trans, atom, atom);

            cxtal->Xcord[at] = atom[0];
            cxtal->Ycord[at] = atom[1];
            cxtal->Zcord[at] = atom[2];
            cxtal->atoms[2*at] = cxtal->atoms[2*asym];
            cxtal->atoms[2*at + 1] = cxtal->atoms[2*asym + 1];
            at++;
        }
    }

}

void apply_all_symmetry_ops(crystal *xtal,
                            molecule *mol,
                            float* mol_Xfrac,
                            float *mol_Yfrac,
                            float *mol_Zfrac,
                            int N,
                            int hall_number)
{
    //get symmetry operations from spglib database
    double translations[192][3];
    int rotations[192][3][3];
    int num_of_operations = spg_get_symmetry_from_database(rotations,
    translations, hall_number);
    int Z = num_of_operations;
    xtal->Z = Z;

    //get lattice_vector inverse and transpose
    float inverse_lattice_vectors[3][3];
    inverse_mat3b3(inverse_lattice_vectors, xtal->lattice_vectors);
    float lattice_vectors_transpose[3][3];
    copy_mat3b3_mat3b3(lattice_vectors_transpose, xtal->lattice_vectors);
    mat3b3_transpose(lattice_vectors_transpose, lattice_vectors_transpose);

    //i loops over all the operations
    for(int i = 0; i < num_of_operations; i++)
    {
        int molecule_index = i*N;

        //rot and trans are ith symmetry operation
        int rot[3][3];
        copy_intmat3b3_intmat3b3bN(rot, rotations, i);
        float trans[3];
        trans[0] = translations[i][0];
        trans[1] = translations[i][1];
        trans[2] = translations[i][2];

        //j loops over all atoms in a molecule
        for(int j = 0; j < N; j++)
        {
            float atomj_array[3];
            atomj_array[0] = mol_Xfrac[j];
            atomj_array[1] = mol_Yfrac[j];
            atomj_array[2] = mol_Zfrac[j];

            //apply operation
            vector3_intmat3b3_multiply(rot, atomj_array, atomj_array);
            vector3_add(trans, atomj_array, atomj_array);

            //temp_array has the fractional cord of jth atom in ith mol
            //convert to cartesian
            vector3_mat3b3_multiply(lattice_vectors_transpose,
                atomj_array, atomj_array);

            //copy to the structure
            xtal->Xcord[molecule_index+j] = atomj_array[0];
            xtal->Ycord[molecule_index+j] = atomj_array[1];
            xtal->Zcord[molecule_index+j] = atomj_array[2];
            xtal->atoms[2*(molecule_index+j)] = (*mol).atoms[2*j];
            xtal->atoms[2*(molecule_index+j)+1] = (*mol).atoms[2*j+1];
        }
    }
    bring_all_molecules_to_first_cell(xtal);
    xtal->num_atoms_in_molecule = mol->num_of_atoms;
}



//////// for layer group ///////////////



void apply_all_lg_symmetry_ops(crystal *xtal,
			       molecule *mol,
			       float* mol_Xfrac,
			       float *mol_Yfrac,
			       float *mol_Zfrac,
			       int N,
			       int hall_number)
{
	//get symmetry operations from spglib database
	double translations[192][3];
	int rotations[192][3][3];

	//this is from spglib, need to be changed
	//int num_of_operations = get_spg_symmetry(hall_number,translations,rotations);
	int num_of_operations = get_lg_symmetry(hall_number,translations,rotations); //changed here

	int Z = num_of_operations;
	xtal->Z = Z;

	//get lattice_vector inverse and transpose
	float inverse_lattice_vectors[3][3];
	inverse_mat3b3(inverse_lattice_vectors, xtal->lattice_vectors);
	float lattice_vectors_transpose[3][3];
	copy_mat3b3_mat3b3(lattice_vectors_transpose, xtal->lattice_vectors);
	mat3b3_transpose(lattice_vectors_transpose, lattice_vectors_transpose);

	//i loops over all the operations
	for(int i = 0; i < num_of_operations; i++)
	{
		int molecule_index = i*N;

		//rot and trans are ith symmetry operation
		int rot[3][3];
		copy_intmat3b3_intmat3b3bN(rot, rotations, i);
		float trans[3];
		trans[0] = translations[i][0];
		trans[1] = translations[i][1];
		trans[2] = translations[i][2];

		//j loops over all atoms in a molecule
		for(int j = 0; j < N; j++)
		{
			float atomj_array[3];
			atomj_array[0] = mol_Xfrac[j];
			atomj_array[1] = mol_Yfrac[j];
			atomj_array[2] = mol_Zfrac[j];

			//apply operation
			vector3_intmat3b3_multiply(rot, atomj_array, atomj_array);
			vector3_add(trans, atomj_array, atomj_array);

			//temp_array has the fractional cord of jth atom in ith mol
			//convert to cartesian
			vector3_mat3b3_multiply(lattice_vectors_transpose,
				atomj_array, atomj_array);

			//copy to the structure
			xtal->Xcord[molecule_index+j] = atomj_array[0];
			xtal->Ycord[molecule_index+j] = atomj_array[1];
			xtal->Zcord[molecule_index+j] = atomj_array[2];
			xtal->atoms[2*(molecule_index+j)] = (*mol).atoms[2*j];
            xtal->atoms[2*(molecule_index+j)+1] = (*mol).atoms[2*j+1];
		}
	}
	//printf("before bring all mol to first Z is : %d\n",(*xtal).Z);
	//fflush(stdout);
	bring_all_molecules_to_first_cell(xtal);
	//printf("after bring all mol to first Z is : %d\n",(*xtal).Z);
	xtal->num_atoms_in_molecule = mol->num_of_atoms;
}



/*function for placing molecule at any position for any spacegroup
 * molecule is automatically aligned according to the requirement
 * of the wyckoff position. for general positions, molecule is rotated
 * randomly. generalizes all the functions written before.
 * assumes symmetry directions to be along cartesian axes of molecule
 * #TODO: create a molecule symm structure to store molecule symmetry
 */
int auto_align_and_generate_at_position(crystal *Xtal,
                                molecule *mol,
                                int hall_number,
                                int spg,
                                int pos_index,
                                COMPATIBLE_SPG compatible_spg)
{

    //declare variables
    //num of atoms in the molecule is N
    int N = mol->num_of_atoms;
    //int Z = Xtal->Z;
    int wyckoff_pos = compatible_spg.allowed_pos[pos_index];
    int order = spg_positions[spg-1].multiplicity[0]/\
        spg_positions[spg-1].multiplicity[wyckoff_pos];

    //get overlap list from the data compatible_spg structure
    int *overlap_list = compatible_spg.pos_overlap_list[pos_index];
    //int num_atoms_in_cell = N*Z;
    //int general_position = 0;
    //From a database, get first coordinates where mol is to be placed
    //#TODO: get rot_mat and trans_vec from a databse of first position.
    float rot_mat[3][3] = {{0,0,1},{0,1,0},{0,0,1}};
    copy_mat3b3_intmat3b3bN(rot_mat,
                            spg_positions[spg-1].first_position_rot,
                            wyckoff_pos);
    float trans_vec[3] = {0, 0 ,0};
    copy_vector3_vector3bN(trans_vec,
                          spg_positions[spg-1].first_position_trans,
                          wyckoff_pos);

    //if general position, rotate molecule randomly
    if ( get_degrees_of_freedom(spg, wyckoff_pos) == 2)
    {
        //create rotation matix around the specified axis with angle psi
        float random_rotation_matrix[3][3];
        generate_random_rotation_matrix(random_rotation_matrix);
        //randomly rotate
        molecule_rotate(mol, random_rotation_matrix);
	get_euler_from_rotation_matrix(random_rotation_matrix, Xtal->euler_angles);
    }
    // place the first molecule
    float inverse_lattice_vectors[3][3], lat_vec_trans[3][3];
    inverse_mat3b3(inverse_lattice_vectors, Xtal->lattice_vectors);
    mat3b3_transpose(inverse_lattice_vectors,inverse_lattice_vectors);
    mat3b3_transpose(lat_vec_trans, Xtal->lattice_vectors);
    float mol_Xfrac[N]; //stores fractional coordinates of first mol
    float mol_Yfrac[N];
    float mol_Zfrac[N];

    float rand_frac_array[3] = {uniform_dist_01(),\
                                uniform_dist_01(),\
                                uniform_dist_01()};
    vector3_mat3b3_multiply(rot_mat,rand_frac_array,rand_frac_array);
    vector3_add(trans_vec,rand_frac_array,rand_frac_array);
    //rand_frac_array is the postion of the first mol

    vector3_mat3b3_multiply(lat_vec_trans, rand_frac_array, Xtal->com_positions);
    
    //compute fractional coordiates of the first position and
    //move molecule to the position
    for(int i = 0; i < N; i++)
    {
        float atom_i[3] = {(*mol).X[i], (*mol).Y[i], (*mol).Z[i]};
        vector3_mat3b3_multiply(inverse_lattice_vectors,
                                atom_i,
                                atom_i);
        vector3_add(atom_i, rand_frac_array, atom_i); //translate
        mol_Xfrac[i] = atom_i[0];
        mol_Yfrac[i] = atom_i[1];
        mol_Zfrac[i] = atom_i[2];

    }
    //now mol_frac has the first mol in fractional coordinates
    // at the given wyckoff position.
    //now apply all the symmetry operations of spg
    apply_all_symmetry_ops(Xtal,
                           mol,
                           mol_Xfrac,
                           mol_Yfrac,
                           mol_Zfrac,
                           N,
                           hall_number);

    //check if inv center or genral position
    int dof = get_degrees_of_freedom(spg, wyckoff_pos);
    if ( dof == 2)
    {
        //if inv centre remove overlap molecules
        if( strcmp(spg_positions[spg-1].site_symmetry[wyckoff_pos], "-1") == 0 )
            combine_close_molecules(Xtal);
        return 1;
    }

    Xtal->spg = spg;
    bring_molecules_to_origin(Xtal);

    int result = align_using_std_orientations(Xtal, mol, hall_number,
        rand_frac_array, overlap_list,
        order, dof, compatible_spg.compatible_axes[wyckoff_pos]);


    if (!result)
    {
        return 0;
    }

    return 1;
}


///////////////// for layer group /////////////////////////

int lg_auto_align_and_generate_at_position(crystal *Xtal,
					   molecule *mol,
					   int hall_number,
					   int spg,
					   int pos_index,
					   COMPATIBLE_SPG compatible_spg)
{

	//declare variables
	//num of atoms in the molecule is N
	int N = mol->num_of_atoms;
	//int Z = Xtal->Z;
	int wyckoff_pos = compatible_spg.allowed_pos[pos_index];  //index of wykoff_pos
	int order = lg_positions[spg-1].multiplicity[0]/\
		lg_positions[spg-1].multiplicity[wyckoff_pos];

	//get overlap list from the data compatible_spg structure
	int *overlap_list = compatible_spg.pos_overlap_list[pos_index];
	//int num_atoms_in_cell = N*Z;
	//int general_position = 0;
	//From a database, get first coordinates where mol is to be placed
	//#TODO: get rot_mat and trans_vec from a databse of first position.
	float rot_mat[3][3] = {{0,0,1},{0,1,0},{0,0,1}};
	copy_mat3b3_intmat3b3bN(rot_mat,lg_positions[spg-1].first_position_rot,wyckoff_pos);
	float trans_vec[3] = {0, 0 ,0};
	copy_vector3_vector3bN(trans_vec,lg_positions[spg-1].first_position_trans,wyckoff_pos);

	//if general position, rotate molecule randomly
	//this is to check site symmetry
	if ( lg_get_degrees_of_freedom(spg, wyckoff_pos) == 2)
	{
		//create rotation matix around the specified axis with angle psi
		float random_rotation_matrix[3][3];
		generate_random_rotation_matrix(random_rotation_matrix);
		//randomly rotate
		molecule_rotate(mol, random_rotation_matrix);
	}


	// place the first molecule
	float inverse_lattice_vectors[3][3];
	inverse_mat3b3(inverse_lattice_vectors, Xtal->lattice_vectors);
	mat3b3_transpose(inverse_lattice_vectors,inverse_lattice_vectors);
	float mol_Xfrac[N]; //stores fractional coordinates of first mol
	float mol_Yfrac[N];
	float mol_Zfrac[N];

	float rand_frac_array[3] = {uniform_dist_01(),\
							uniform_dist_01(),\
								uniform_dist_01()};
	vector3_mat3b3_multiply(rot_mat,rand_frac_array,rand_frac_array);
	vector3_add(trans_vec,rand_frac_array,rand_frac_array);
	//rand_frac_array is the postion of the first mol

	//compute fractional coordiates of the first position and
	//move molecule to the position
	for(int i = 0; i < N; i++)
	{
		float atom_i[3] = {(*mol).X[i], (*mol).Y[i], (*mol).Z[i]};
		vector3_mat3b3_multiply(inverse_lattice_vectors,
								atom_i,
								atom_i);
		vector3_add(atom_i, rand_frac_array, atom_i); //translate
		mol_Xfrac[i] = atom_i[0];
		mol_Yfrac[i] = atom_i[1];
		mol_Zfrac[i] = atom_i[2];

	}
	//now mol_frac has the first mol in fractional coordinates
	// at the given wyckoff position.
	//printf("Z after align first molecule: %d\n",(*Xtal).Z);
	//fflush(stdout);

	//now apply all the symmetry operations of spg
	apply_all_lg_symmetry_ops(Xtal,
				  mol,
			          mol_Xfrac,
			          mol_Yfrac,
				  mol_Zfrac,
				  N,
				  hall_number);

	//check if inv center or genral position
	//printf("Z after apply_all_symmetry_ops: %d\n",(*Xtal).Z);
	//fflush(stdout);

	//this is to check site symmetry
	int dof = lg_get_degrees_of_freedom(spg, wyckoff_pos);
	if ( dof == 2)
	{
		//if inv centre remove overlap molecules
		if( strcmp(lg_positions[spg-1].site_symmetry[wyckoff_pos], "-1") == 0 )
			combine_close_molecules(Xtal);

		//printf("dof is 2,Z after apply_all_symmetry_ops: %d\n",(*Xtal).Z);
		//fflush(stdout);
		return 1;
	}

	Xtal->spg = spg;
	bring_molecules_to_origin(Xtal);

	int result = lg_align_using_std_orientations(Xtal, mol, hall_number,
		rand_frac_array, overlap_list,
		order, dof, compatible_spg.compatible_axes[wyckoff_pos]);

	/*
	printf("Z at end of auto_align_and_generate_at_position: %d\n",(*Xtal).Z);
	fflush(stdout);
	*/

	if (!result)
	{
		//printf("spg: %d , pos: %d unable to align\n", spg, wyckoff_pos);
		return 0;
	}

    //printf("spg: %d , pos: %d aligned\n", spg, wyckoff_pos);
	return 1;

}



int align_using_std_orientations(crystal* xtal_1,
                                molecule* mol,
                                int hall_number,
                                float first_com[3],
                                int overlap_list[],
                                int len_overlap_list,
                                int dof,
                                COMPATIBLE_AXES compatible_axes)
{
    //take two pairs
    //rotate molecule to average position
    int N = mol->num_of_atoms;
    int Z = xtal_1->Z;
    float rotation_matrix[3][3];
    float mol_Xfrac[N]; //stores fractional ...
    float mol_Yfrac[N]; //coordinates of first mol
    float mol_Zfrac[N];
    float lattice_vectors[3][3];
    float inv_lattice_vectors[3][3];
    copy_mat3b3_mat3b3(lattice_vectors, xtal_1->lattice_vectors);
    inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
    mat3b3_transpose(inv_lattice_vectors, inv_lattice_vectors);

    float *mol_axes = compatible_axes.usable_mol_axes;
    float *viewing_directions = compatible_axes.usable_view_dir;
    int comb = compatible_axes.num_combinations;

    //create a temp xtal
    crystal *xtal = (crystal *)malloc(sizeof(crystal) );
    allocate_xtal(xtal, Z, N);
    copy_xtal(xtal, xtal_1);

    //rotate about axis if allowed
    float random_rotation_about_axis[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    int i = rand_r(seed2) % comb;

    float mol_axis[3] = {*(mol_axes+3*i + 0),
                         *(mol_axes+3*i + 1),
                         *(mol_axes+3*i + 2) };
    float view_dir[3] = { *(viewing_directions + 3*i + 0),
                          *(viewing_directions + 3*i + 1),
                          *(viewing_directions + 3*i + 2)};
    normalise_vector3(view_dir);
    normalise_vector3(mol_axis);

    //create rotation matrix
    rotation_matrix_from_vectors(rotation_matrix, mol_axis, view_dir);
    //rotate and store first molecule in fractional coord
    //in mol_frac
    if(dof == 1)
    {
        float psi = 2*PI*uniform_dist_01();
        rotation_mat_around_axis(random_rotation_about_axis, view_dir, psi);
    }

    if(dof == 0)
    {
        float axis[3];
        int draw = rand_r(seed2) % 3;

        if (draw == 1)
        {
            float x[3] = {1, 0, 0};
            copy_vector3_vector3(axis, x);
        }
        else if (draw == 2)
        {
            float y[3] = {0, 1, 0};
            copy_vector3_vector3(axis, y);
        }
        else
        {
            float z[3] = {0, 0, 1};
            copy_vector3_vector3(axis, z);
        }
        rotation_mat_around_axis(random_rotation_about_axis, axis, PI/2);
    }

    for(int z = 0; z < N; z++)
    {
        float temp[3] = {xtal_1->Xcord[z] ,
                         xtal_1->Ycord[z] ,
                         xtal_1->Zcord[z] };
        vector3_mat3b3_multiply(rotation_matrix, temp, temp);

        vector3_mat3b3_multiply(random_rotation_about_axis, temp, temp);
        //vector3_add(temp,com,temp);
        //convert to frac
        vector3_mat3b3_multiply(inv_lattice_vectors,
                                temp,
                                temp );
        vector3_add(temp, first_com, temp);
        mol_Xfrac[z] = temp[0];
        mol_Yfrac[z] = temp[1];
        mol_Zfrac[z] = temp[2];
    }

    apply_all_symmetry_ops( xtal,
                            mol,
                            mol_Xfrac,
                            mol_Yfrac,
                            mol_Zfrac,
                            N,
                            hall_number);

    //bring_molecules_to_origin(xtal);

    int result = check_overlap_xtal_cartesian(xtal,
                                        overlap_list,
                                        len_overlap_list,
                                        N);

    if (result)
    {
        int Z_gen = xtal->Z;
        combine_close_molecules(xtal);
        int Z_return = xtal->Z;
        if(Z_return != Z_gen/len_overlap_list)
            {goto failed;}

        copy_xtal(xtal_1, xtal);
        goto success;
    }

    failed:
        free_xtal(xtal);
        return 0;
    success:
        free_xtal(xtal);
        return 1;

}


/////////////  for layer group ////////////////////


int lg_align_using_std_orientations(crystal* xtal_1,
				    molecule* mol,
				    int hall_number,
				    float first_com[3],
				    int overlap_list[],
				    int len_overlap_list,
			            int dof,
				    COMPATIBLE_AXES compatible_axes)
{
	//take two pairs
	//rotate molecule to average position
	int N = mol->num_of_atoms;
	int Z = xtal_1->Z;
	float rotation_matrix[3][3];
	float mol_Xfrac[N]; //stores fractional ...
	float mol_Yfrac[N]; //coordinates of first mol
	float mol_Zfrac[N];
	float lattice_vectors[3][3];
	float inv_lattice_vectors[3][3];
	copy_mat3b3_mat3b3(lattice_vectors, xtal_1->lattice_vectors);
	inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
	mat3b3_transpose(inv_lattice_vectors, inv_lattice_vectors);

	float *mol_axes = compatible_axes.usable_mol_axes;
	float *viewing_directions = compatible_axes.usable_view_dir;
	int comb = compatible_axes.num_combinations;

	//create a temp xtal
	crystal *xtal = (crystal *)malloc(sizeof(crystal) );
	allocate_xtal(xtal, Z, N);
	copy_xtal(xtal, xtal_1);
	//debug
	//int spg = xtal_1->spg;
	//generate_lattice(xtal->lattice_vectors, spg, 60, 120, 1000000);
	//generate_fake_lattice(xtal->lattice_vectors, spg);

	//rotate about axis if allowed
	float random_rotation_about_axis[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

	int i = rand_r(seed2) % comb;

	float mol_axis[3] = {*(mol_axes+3*i + 0),
						 *(mol_axes+3*i + 1),
						 *(mol_axes+3*i + 2) };
	float view_dir[3] = { *(viewing_directions + 3*i + 0),
						  *(viewing_directions + 3*i + 1),
						  *(viewing_directions + 3*i + 2)};
	normalise_vector3(view_dir);
	normalise_vector3(mol_axis);

	//create rotation matrix
	rotation_matrix_from_vectors(rotation_matrix, mol_axis, view_dir);
	//rotate and store first molecule in fractional coord
	//in mol_frac
	if(dof == 1)
	{
		float psi = 2*PI*uniform_dist_01();
		rotation_mat_around_axis(random_rotation_about_axis, view_dir, psi);
	}

	if(dof == 0)
	{
		float axis[3];
		int draw = rand_r(seed2) % 3;

		if (draw == 1)
		{
			float x[3] = {1, 0, 0};
			copy_vector3_vector3(axis, x);
		}
		else if (draw == 2)
		{
			float y[3] = {0, 1, 0};
			copy_vector3_vector3(axis, y);
		}
		else
		{
			float z[3] = {0, 0, 1};
			copy_vector3_vector3(axis, z);
		}
		rotation_mat_around_axis(random_rotation_about_axis, axis, PI/2);
	}

	for(int z = 0; z < N; z++)
	{
		float temp[3] = {xtal_1->Xcord[z] ,
						 xtal_1->Ycord[z] ,
						 xtal_1->Zcord[z] };
		vector3_mat3b3_multiply(rotation_matrix, temp, temp);

		vector3_mat3b3_multiply(random_rotation_about_axis, temp, temp);
		//vector3_add(temp,com,temp);
		//convert to frac
		vector3_mat3b3_multiply(inv_lattice_vectors,
								temp,
								temp );
		vector3_add(temp, first_com, temp);
		mol_Xfrac[z] = temp[0];
		mol_Yfrac[z] = temp[1];
		mol_Zfrac[z] = temp[2];
	}




	apply_all_lg_symmetry_ops(xtal,
				  mol,
				  mol_Xfrac,
				  mol_Yfrac,
				  mol_Zfrac,
				  N,
				  hall_number);

	//bring_molecules_to_origin(xtal);




	int result = check_overlap_xtal_cartesian(xtal,
						  overlap_list,
						  len_overlap_list,
						  N);

	if (result)
	{
		int Z_gen = xtal->Z;
		//printf("i=%d, j1=%d, j2=%d, k1=%d, k2=%d \n", i,j1,j2,k1,k2);
		//remove_close_molecules(xtal);
		combine_close_molecules(xtal);
		int Z_return = xtal->Z;
		//len_overlap_list == order, Z/len = multiplicity
		if(Z_return != Z_gen/len_overlap_list)
			{goto failed;}

		copy_xtal(xtal_1, xtal);
		goto success;
	}

	failed:
		free_xtal(xtal);
		return 0;
	success:
		free_xtal(xtal);
		return 1;

}


int get_degrees_of_freedom(int spg, int pos)
{
    char pg[6];
    strcpy(pg, spg_positions[spg-1].site_symmetry[pos]);
    //general_position or inversion_centre: complete freedom
    if(!strcmp(pg,"-1") || !strcmp(pg,"1") )
        return 2;
    else if (!strcmp(pg,"2") || !strcmp(pg,"3") || !strcmp(pg,"m") ||
             !strcmp(pg,"4") || !strcmp(pg,"6") || !strcmp(pg, "2/m") ||
             !strcmp(pg,"-4")|| !strcmp(pg,"4/m") || !strcmp(pg,"6/m") ||
             !strcmp(pg,"-3") || !strcmp(pg,"-6")
             )
        {
            return 1;
        }
    else
        return 0;

}


/////////////////////  for layer group //////////////////////


int lg_get_degrees_of_freedom(int spg, int pos)
{
    char pg[6];
    strcpy(pg, lg_positions[spg-1].site_symmetry[pos]);
    //general_position or inversion_centre: complete freedom
    if(!strcmp(pg,"-1") || !strcmp(pg,"1") )
        return 2;
    else if (!strcmp(pg,"2") || !strcmp(pg,"3") || !strcmp(pg,"m") ||
             !strcmp(pg,"4") || !strcmp(pg,"6") || !strcmp(pg, "2/m") ||
             !strcmp(pg,"-4")|| !strcmp(pg,"4/m") || !strcmp(pg,"6/m") ||
             !strcmp(pg,"-3") || !strcmp(pg,"-6")
             )
        {
            return 1;
        }
    else
        return 0;

}
