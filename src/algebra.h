#ifndef ALGEBRA_H
#define ALGEBRA_H

void vector3_add(float a[3], float b[3], float c[3]);
void vector3_mat3b3_multiply(float a[3][3], float b[3], float c[3]);
void mat3b3_mat3b3_multiply(float a[3][3], float b[3][3], float c[3][3]);
void vector3_intmat3b3_multiply(int a[3][3], float b[3], float c[3]);
void rotation_mat_around_axis(float rot[3][3], float axis[3], float psi);
void print_mat3b3(float mat[3][3]);
void print_vec3(float vec[3]);
float det_mat3b3(float a[3][3]);
void copy_vector3_vector3(float a[3], float b[3]);
void copy_mat3b3_mat3b3(float a[3][3], float b[3][3]);
void copy_mat3b3_mat3b3bN(float a[3][3], float b[][3][3], int index);
void copy_mat3b3bN_mat3b3(float b[][3][3], float a[3][3], int index);
void copy_intmat3b3_intmat3b3bN(int a[3][3], int b[][3][3], int index);
void copy_mat3b3_intmat3b3bN(float a[3][3], const int b[][3][3], int index);
void copy_vector3_vector3bN(float a[3], const float b[][3], int index);
void inverse_mat3b3(float Tinv[3][3], float T[3][3]);
int hall_number_from_spg(int spg);
void copy_floatmat3b3_intmat3b3(float a[3][3], int b[3][3]);
void mat3b3_transpose(float b_trans[3][3], float b[3][3]);
void print_mat3b3bN(float a[][3][3], int N);
void copy_vector3_mat3b3(float a[3], float b[3][3],int index);
void vector3_inverse(float a[3]);
float vector3_norm(float a[3]);
void vector3_subtract(float a[3], float b[3], float c[3]);
void cross_vector3_vector3(float cross[3], float a[3], float b[3]);
void normalise_vector3(float a[3]);
float dot_vector3_vector3(float a[3], float b[3]);
int check_vec3_isNull(float a[3], float tol);
float cart_dist(float p1[3], float p2[3]);
void array_shuffler_1(float *a, int len );
void array_shuffler_2(float a[][3], int len );
void generate_random_rotation_matrix( float rotation_matrix[3][3] );
void generate_random_translation_vector(float trans[3]);
int are_equal_floats(float a, float b, float ftol);
void vector3_frac(float a[3]);

void copy_intmat3b3bN_intmat3b3(int b[][3][3], int a[3][3], int index);//added here
void copy_intmat3b3_constintmat3b3bN(int a[3][3], int const b[][3][3], int index);//added here
void copy_vector3bN_vector3(float a[3], float b[][3], int index); //added here
void copy_doubvector3_vector3bN( double a[3], const double b[][3], int index);//added here
void copy_doubvector3bN_vector3(double a[3], double b[][3], int index); //added here
int get_lg_symmetry(int hall_number,double translations [192] [3],int rotations[192][3][3]);//added here

#endif
