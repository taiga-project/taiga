// velocity Verlet

#include "solvers.cuh"

__device__ void solve_diffeq_by_verlet(double *X, double *a, double *B, double *E, double eperm, double timestep){
    double R = X[0] + X[3]*timestep + 0.5*a[0]*timestep*timestep;
    double Z = X[1] + X[4]*timestep + 0.5*a[1]*timestep*timestep;
    double T = X[2] + X[5]*timestep + 0.5*a[2]*timestep*timestep;
    double new_a[3];
    (*get_acceleration_from_lorentz_force)(&X[3], new_a, B, E, eperm);
    X[0] = R;
    X[1] = Z;
    X[2] = T;
    X[3] += (a[0] + new_a[0])*timestep*0.5;
    X[4] += (a[1] + new_a[1])*timestep*0.5;
    X[5] += (a[2] + new_a[2])*timestep*0.5;
    memcpy(a, new_a, 3*sizeof(double));
}
