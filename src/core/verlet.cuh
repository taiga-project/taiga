#ifndef VERLET_CUH
#define VERLET_CUH

__device__ void calculate_verlet_x(double *X, double *B, double *E,double eperm, double timestep);
__device__ void calculate_verlet_v(double *X, double *B, double *E, double *E_prev, double eperm, double dt_per_2);
__device__ void solve_diffeq_by_verlet(double *X, double eperm, double timestep,
                                       TaigaCommons *c, bool is_electric_field_on,
                                       int *local_spline_indices,
                                       double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                       double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                       double *local_psi_n);
#endif //VERLET_CUH
