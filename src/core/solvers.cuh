#ifndef SOLVERS_CUH
#define SOLVERS_CUH

__device__ void (*get_acceleration_from_lorentz_force)(double *a, double *v, double *B, double *E, double eperm);

#endif
