#ifndef RUNGE_KUTTA_NYSTROM_CUH
#define RUNGE_KUTTA_NYSTROM_CUH

__device__ double solve_diffeq_by_rkn(double *X, double eperm, double timestep,
                                      TaigaCommons *c, bool is_electric_field_on,
                                      int *local_spline_indices,
                                      double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                      double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                      double *local_spline_psi_n);

#endif //RUNGE_KUTTA_NYSTROM_CUH
