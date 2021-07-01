// velocity Verlet
#include "solvers.cuh"

__device__ double cross(double *u, double *v, int index) {
    switch(index){
        case 0:
            return u[1]*v[2] - u[2]*v[1];
        case 1:
            return u[2]*v[0] - u[0]*v[2];
        case 2:
            return u[0]*v[1] - u[1]*v[0];
    }
}

__device__ void solve_diffeq_by_verlet(double *X, double *X_prev, double *B, double *E,  double *E_prev, double eperm, double timestep) {
    int i;
    double a[3];

    double t[3];
    double v_minus[3];
    double v_star[3];
    double v_plus[3];

    double t_sq = 0;
    double s__t;

    double dt__2 = 0.5 * timestep;

    for (i = 0; i < 3; ++i) {
        t[i] = eperm * dt__2 * B[i];
        t_sq += t[i] * t[i];
    }
    s__t = 2.0/(1.0+t_sq);

    (*get_acceleration_from_lorentz_force)(a, &X_prev[3], B, E, eperm);

    for (i = 0; i < 3; ++i) {
        X[i] += timestep * (X_prev[i+3] + dt__2*a[i]);
        E[i] = 0.5 * (E[i] + E_prev[i]);
        v_minus[i] = X_prev[i+3] + dt__2*E[i];
    }

    for (i = 0; i < 3; ++i) {
        v_star[i] = v_minus[i] + cross(v_minus, t, i);
    }

    for (i = 0; i < 3; ++i) {
        v_plus[i] = v_minus[i] + s__t * cross(v_star, t, i);
        X[i+3] = v_plus[i] + dt__2 * E[i];
    }

}
