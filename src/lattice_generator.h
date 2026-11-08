#ifndef LATTICE_GENERATOR_H
#define LATTICE_GENERATOR_H

/* A failed principal-component draw (infeasible volume or exhausted bounded
 * rejection sampling) returns nonfinite diagonal entries. Callers must reject
 * that candidate before coordinate conversion. Public void signatures remain.
 */
void gen_triclinic_lattice(float lattice_vector[3][3],
                           float target_volume,
                           float norm_std,
                           float angle_std);

void gen_monoclinic_lattice(float lattice_vector[3][3],
                            float target_volume,
                            float norm_std,
                            float angle_std);

void gen_orthorhombic_lattice(float lattice_vector[3][3],
                              float target_volume,
                              float norm_std);

void gen_tetragonal_lattice(float lattice_vector[3][3],
                            float target_volume,
                            float norm_std);

void gen_hexagonal_lattice(float lattice_vector[3][3],
                           float target_volume,
                           float norm_std);

void gen_cubic_lattice(float lattice_vector[3][3],
                       float target_volume);

void generate_lattice(float lattice_vector[3][3],
                      int spg,
                      float norm_std,
                      float angle_std,
                      float target_volume);

void generate_fake_lattice(float lattice_vector[3][3], int spg);

#endif

