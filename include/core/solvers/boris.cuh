#ifndef BORIS_CUH
#define BORIS_CUH

__device__ void calculate_boris_v(double *X, double *B, double *E, double *E_prev, double eperm, double dt_per_2);

#endif //BORIS_CUH
