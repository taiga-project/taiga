#ifndef RK4_CUH
#define RK4_CUH

__device__ void calculate_runge_kutta_coeff(double *X,
                                            double *S, double *S_prev, double rk_weight,
                                            double *B, double *E,
                                            double eperm, double timestep);

__device__ double solve_diffeq_by_rk4(double *X, double eperm, double timestep,
                                      TaigaCommons *c, bool is_electric_field_on,
                                      int *local_spline_indices,
                                      double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                      double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                      double *local_spline_psi_n);

#endif //RK4_CUH
