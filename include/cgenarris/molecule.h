#ifndef _MOLECULE_H
#define _MOLECULE_H

typedef struct
{
    char *atoms;
    float *X;
    float *Y;
    float *Z;
    int num_of_atoms;
}molecule;

#endif