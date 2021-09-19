#include "core/solvers/boris.cuh"

__device__ void calculate_boris_v(double *X, double *B, double *E, double *E_prev, double eperm, double timestep) {
    int i;
    double t[3];
    double v_minus[3];
    double v_star[3];
    double v_plus[3];
    double t_square = 0;
    double s_per_t;
    double dt_per_2 = timestep / 2.0;

    for (i = 0; i < 3; ++i) {
        t[i] = eperm * dt_per_2 * B[i];
        t_square += t[i] * t[i];
    }
    s_per_t = 2.0 / (1.0 + t_square);

    for (i = 0; i < 3; ++i) {
        E[i] = 0.5 * (E[i] + E_prev[i]);
        v_minus[i] = X[i + 3] + eperm * dt_per_2 * E[i];
    }

    for (i = 0; i < 3; ++i) {
        v_star[i] = v_minus[i] + cross(v_minus, t, i);
    }

    for (i = 0; i < 3; ++i) {
        v_plus[i] = v_minus[i] + s_per_t * cross(v_star, t, i);
        X[i + 3] = v_plus[i] + eperm * dt_per_2 * E[i];
    }
}
