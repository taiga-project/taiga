#ifndef VERLET_CUH
#define VERLET_CUH

__device__ void calculate_verlet_x(double *X, double *B, double *E,double eperm, double timestep);
__device__ void calculate_verlet_v(double *X, double *B, double *E, double *E_prev, double eperm, double dt_per_2);
__device__ void solve_diffeq_by_verlet(double *X, double *B, double *E, double *E_prev, double eperm, double timestep);

#endif //VERLET_CUH
