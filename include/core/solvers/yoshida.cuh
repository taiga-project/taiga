#ifndef YOSHIDA_CUH
#define YOSHIDA_CUH

__device__ void calculate_yoshida_x(double c, double *X, double timestep);
__device__ void calculate_yoshida_v(double d, double *X,
                                    double *B, double *E, double *E_prev,
                                    double eperm, double timestep);
__device__ double solve_diffeq_by_yoshida(double *X, double eperm, double timestep,
                                          TaigaCommons *c, bool is_electric_field_on,
                                          int *local_spline_indices,
                                          double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                          double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                          double *local_spline_psi_n);

#endif //YOSHIDA_CUH
