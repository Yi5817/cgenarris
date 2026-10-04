#ifndef __SYMMETRIZATION_H__
#define __SYMMETRIZATION_H__

void symmetrize_state(double *state, int *invert, const int nmol, const int spg);
void symmetrize_gradient(double *grad, int *invert, const int nmol, const int spg);
#endif
