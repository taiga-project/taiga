// Runge--Kutta method

#include "core/solvers/rk4.cuh"
#include "core/solvers/solvers.cuh"
#include "core/localise_field.cuh"

__device__ void calculate_runge_kutta_coeff(double *X,
                                            double *S, double *S_prev, double rk_weight,
                                            double *B, double *E,
                                            double eperm, double timestep){
    int i;
    for (i=0; i<3; ++i){
        S[i] = X[i+3] + rk_weight * S_prev[i+3];
    }
    (*get_acceleration_from_lorentz_force)(&S[3], S, B, E, eperm);
    
    for (i=0; i<6; ++i){
        S[i] *= timestep;
    }
}

__device__ double solve_diffeq_by_rk4(double *X, double eperm, double timestep,
                                    TaigaCommons *c, bool is_electric_field_on,
                                    int *local_spline_indices,
                                    double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                    double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                    double *local_spline_psi_n){
    double B[3], E[3];
    double local_psi_n;
    local_psi_n = get_local_field(X, B, E, c, is_electric_field_on,
                                  local_spline_indices,
                                  local_spline_brad, local_spline_bz, local_spline_btor,
                                  local_spline_erad, local_spline_ez, local_spline_etor,
                                  local_spline_psi_n);

    double S1[6], S2[6], S3[6], S4[6];
    
    calculate_runge_kutta_coeff(X, S1, X,  0.0, B, E, eperm ,timestep);
    calculate_runge_kutta_coeff(X, S2, S1, 0.5, B, E, eperm ,timestep);
    calculate_runge_kutta_coeff(X, S3, S2, 0.5, B, E, eperm ,timestep);
    calculate_runge_kutta_coeff(X, S4, S3, 1.0, B, E, eperm ,timestep);

    int i;
    for(i=0; i<6; ++i){
        X[i] += (S1[i] + 2*S2[i] + 2*S3[i] + S4[i])/6;
    }

    return local_psi_n;
}
