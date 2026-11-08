#ifndef RANDOMGEN_H
#define RANDOMGEN_H

void init_genrand(unsigned int s);
void init_by_array(unsigned long init_key[], int key_length);
void init_genrand_rank(unsigned int s, int rank);
unsigned long genrand_int32(void);
float uniform_dist_01 (void);
float normal_dist_01 (void);
float normal_dist_ab ( float a, float b);

extern unsigned int *seed2;
int get_random_int(void);

#endif
