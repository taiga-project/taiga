// Runge--Kutta method

#include "solvers.cuh"

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

__device__ void solve_diffeq_by_rk4(double *X, double *B, double *E, double *E_prev, double eperm, double timestep){
    double S1[6], S2[6], S3[6], S4[6];
    
    calculate_runge_kutta_coeff(X, S1, X,  0.0, B, E, eperm ,timestep);
    calculate_runge_kutta_coeff(X, S2, S1, 0.5, B, E, eperm ,timestep);
    calculate_runge_kutta_coeff(X, S3, S2, 0.5, B, E, eperm ,timestep);
    calculate_runge_kutta_coeff(X, S4, S3, 1.0, B, E, eperm ,timestep);

    int i;
    for(i=0; i<6; ++i){
        X[i] += (S1[i] + 2*S2[i] + 2*S3[i] + S4[i])/6;
    }
}