// Runge--Kutta--Nystrom method
#include "core/solvers/runge_kutta_nystrom.cuh"
#include "core/solvers/solvers.cuh"
#include "core/localise_field.cuh"

__device__ double solve_diffeq_by_rkn(double *X, double eperm, double timestep,
                                    TaigaCommons *c, bool is_electric_field_on,
                                    int *local_spline_indices,
                                    double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                    double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                    double *local_spline_psi_n){

    double K1[3], K2[3], K3[3], K4[3];
    double B[3], E[3];
    double x[3], v[3];
    double local_psi_n;
    int i;
    double timestep_square = timestep*timestep;

    //K1
    get_local_field(X, B, E, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    (*get_acceleration_from_lorentz_force)(K1, &X[3], B, E, eperm);

    //K2
    for(i=0; i<3; ++i){
        x[i] = X[i] + timestep / 2.0 * X[i+3] + timestep_square / 8.0 * K1[i];
        v[i] = X[i+3] + timestep / 2.0 * K1[i];
    }
    get_local_field(x, B, E, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    (*get_acceleration_from_lorentz_force)(K2, v, B, E, eperm);

    //K3
    for(i=0; i<3; ++i){
        v[i] = X[i+3] + timestep / 2.0 * K2[i];
    }
    get_local_field(x, B, E, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    (*get_acceleration_from_lorentz_force)(K3, v, B, E, eperm);

    //K4
    for(i=0; i<3; ++i){
        x[i] = X[i] + timestep * X[i+3] + timestep_square / 2.0 * K3[i];
        v[i] = X[i+3] + timestep * K3[i];
    }
    local_psi_n = get_local_field(x, B, E, c, is_electric_field_on,
                                  local_spline_indices,
                                  local_spline_brad, local_spline_bz, local_spline_btor,
                                  local_spline_erad, local_spline_ez, local_spline_etor,
                                  local_spline_psi_n);
    (*get_acceleration_from_lorentz_force)(K4, v, B, E, eperm);

    // end
    for(i=0; i<3; ++i){
        X[i] += timestep * X[i+3] + (K1[i] + K2[i] + K3[i]) * timestep * timestep / 6.0;
        X[i+3] += (K1[i] + 2 * K2[i] + 2 * K3[i] + K4[i]) * timestep / 6.0;
    }
    return local_psi_n;
}
