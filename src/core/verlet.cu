// velocity Verlet with Boris algorithm
#include "verlet.cuh"
#include "solvers.cuh"
#include "maths.cuh"

__device__ void calculate_verlet_x(double *X, double *B, double *E,double eperm, double timestep) {
    int i;
    double a[3];
    double dt_per_2 = 0.5 * timestep;

    (*get_acceleration_from_lorentz_force)(a, &X[3], B, E, eperm);

    for (i = 0; i < 3; ++i) {
        X[i] += timestep * (X[i + 3] + dt_per_2 * a[i]);
    }
}

__device__ void calculate_verlet_v(double *X, double *B, double *E, double *E_prev, double eperm, double dt_per_2) {
    int i;
    double t[3];
    double v_minus[3];
    double v_star[3];
    double v_plus[3];
    double t_square = 0;
    double s_per_t;

    for (i = 0; i < 3; ++i) {
        t[i] = eperm * dt_per_2 * B[i];
        t_square += t[i] * t[i];
    }
    s_per_t = 2.0 / (1.0 + t_square);

    for (i = 0; i < 3; ++i) {
        E[i] = 0.5 * (E[i] + E_prev[i]);
        v_minus[i] = X[i+3] + dt_per_2 * E[i];
    }

    for (i = 0; i < 3; ++i) {
        v_star[i] = v_minus[i] + cross(v_minus, t, i);
    }

    for (i = 0; i < 3; ++i) {
        v_plus[i] = v_minus[i] + s_per_t * cross(v_star, t, i);
        X[i+3] = v_plus[i] + dt_per_2 * E[i];
    }
}

__device__ double solve_diffeq_by_verlet(double *X, double eperm, double timestep,
                                       TaigaCommons *c, bool is_electric_field_on,
                                       int *local_spline_indices,
                                       double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                       double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                       double *local_spline_psi_n){
    double B[3], E[3], E_prev[3];
    double dt_per_2;
    double local_psi_n;

    get_local_field(X, B, E_prev, c, is_electric_field_on,
                    local_spline_indices,
                    local_spline_brad, local_spline_bz, local_spline_btor,
                    local_spline_erad, local_spline_ez, local_spline_etor,
                    local_spline_psi_n);
    calculate_verlet_x(X, B, E_prev, eperm, timestep);
    local_psi_n = get_local_field(X, B, E, c, is_electric_field_on,
                                  local_spline_indices,
                                  local_spline_brad, local_spline_bz, local_spline_btor,
                                  local_spline_erad, local_spline_ez, local_spline_etor,
                                  local_spline_psi_n);
    dt_per_2 = 0.5 * timestep;
    calculate_verlet_v(X, B, E, E_prev, eperm, dt_per_2);
    return local_psi_n;
}