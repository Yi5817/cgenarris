#ifndef LATTICE_GENERATOR_LAYER_H
#define LATTICE_GENERATOR_LAYER_H

//void gen_triclinic_lattice(float lattice_vector[3][3], 
//	float target_volume, float max_angle, float min_angle);
//void gen_monoclinic_lattice(float lattice_vector[3][3], 
//	float target_volume, float max_angle, float min_angle);
/*
void gen_orthorhombic_lattice(float lattice_vector[3][3], 
	float target_volume);
void gen_tetragonal_lattice(float lattice_vector[3][3],
	float target_volume);
void gen_hexagonal_lattice(float lattice_vector[3][3], 
	float target_volume);
void gen_cubic_lattice(float lattice_vector[3][3],
	float target_volume);

*/
void generate_oblique_triclinic(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_oblique_monoclinic(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_rectangle_monoclinic(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_rectangle_orthorhombic(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
                                                        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_hexagonal_hexagonal(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
                                                       float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_square_tetragonal(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
                                                        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);


void generate_layer_lattice(int *all_substrate_combo,float lattice_vector[3][3], int spg,
 float max_angle, float min_angle, float target_volume, float lattice_vector_2d [2][3],int num_combo,
 float interface_area_mean,float interface_area_std,int volume_multiplier,int SET_INTERFACE_AREA);

void generate_fake_layer_lattice(float lattice_vector[3][3], int spg);

void standardise_layer_lattice(float lattice[3][3], int spg);

int generate_substrate_lattice_combs(int *all_substrate_combo, float lattice_vector_2d[2][3],
				 float target_volume, float max_angle,float min_angle );


void get_lower_triangle(float lattice_vector[3][3],float lower_lattice_vector[3][3]);
#endif

