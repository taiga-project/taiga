// Yoshida integrator

#include "core/solvers/boris.cuh"

__device__ void calculate_yoshida_x(double c, double *X, double *B, double *E,double eperm, double timestep) {
    int i;
    double c_dt =  c * timestep;
    for (i = 0; i < 3; ++i) {
        X[i] += c_dt * X[i+3];
    }
}

__device__ void calculate_yoshida_v(double d, double *X,
                                    double *B, double *E, double *E_prev,
                                    double eperm, double timestep){
    calculate_boris_v(X, B, E, E_prev, eperm, d/2.0 * timestep);
}

__device__ double solve_diffeq_by_yoshida(double *X, double eperm, double timestep,
                                        TaigaCommons *c, bool is_electric_field_on,
                                        int *local_spline_indices,
                                        double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                        double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                        double *local_spline_psi_n){
    double B[3], E[3], E_prev[3];
    double local_psi_n;
    double cbrt2 = cbrt(2.0);
    double w1 = 1.0 / (2.0 - cbrt2);
    double w0 = - cbrt2 * w1;
    double c1 = 0.5 * w1;
    double c2 = 0.5 * (w0 + w1);

    get_local_field(X, B, E_prev, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    calculate_yoshida_x(c1, X, B, E, eperm, timestep);
    get_local_field(X, B, E, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    calculate_yoshida_v(w1, X, B, E, E_prev, eperm, timestep);
    get_local_field(X, B, E_prev, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    calculate_yoshida_x(c2, X, B, E, eperm, timestep);
    get_local_field(X, B, E, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    calculate_yoshida_v(w0, X, B, E, E_prev, eperm, timestep);
    get_local_field(X, B, E_prev, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    calculate_yoshida_x(c2, X, B, E, eperm, timestep);
    get_local_field(X, B, E, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    calculate_yoshida_v(w1, X, B, E, E_prev, eperm, timestep);
    local_psi_n = get_local_field(X, B, E_prev, c, is_electric_field_on,
                                  local_spline_indices,
                                  local_spline_brad, local_spline_bz, local_spline_btor,
                                  local_spline_erad, local_spline_ez, local_spline_etor,
                                  local_spline_psi_n);
    calculate_yoshida_x(c1, X, B, E, eperm, timestep);
    return local_psi_n;
}
