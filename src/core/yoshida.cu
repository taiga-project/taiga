// Yoshida integrator

#include "solvers.cuh"
#include "verlet.cuh"

__device__ void calculate_yoshida_x(double c, double *X, double timestep){
    double c_dt = c * timestep;
    int i;
    for (i = 0; i < 3; ++i) {
        X[i] += c_dt * X[i + 3];
    }
}

__device__ void calculate_yoshida_v(double d, double *X,
                                    double *B, double *E, double *E_prev,
                                    double eperm, double timestep){
    calculate_verlet_v(X, B, E, E_prev, eperm, d * timestep);
}

__device__ void solve_diffeq_by_yoshida(double *X, double *B, double *E, double *E_prev, double eperm, double timestep){
    double cbrt2 = cbrt(2.0);
    double w1 = 1.0 / (2.0 - cbrt2);
    double w0 = - cbrt2 * w1;
    double c1 = 0.5 * w1;
    double c2 = 0.5 * (w0 + w1);
    calculate_yoshida_x(c1, X, timestep);
    calculate_yoshida_v(w1, X, B, E, E_prev, eperm, timestep);
    calculate_yoshida_x(c2, X, timestep);
    calculate_yoshida_v(w0, X, B, E, E_prev, eperm, timestep);
    calculate_yoshida_x(c2, X, timestep);
    calculate_yoshida_v(w1, X, B, E, E_prev, eperm, timestep);
    calculate_yoshida_x(c1, X, timestep);
}
