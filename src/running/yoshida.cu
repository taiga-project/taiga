// Yoshida integrator
// velocity Verlet

#include "lorentz.cu"

__device__ void calculate_yoshida_x(double c, double *X_new, double *X, double *a, double timestep){
    double c_dt = c*timestep;
    X_new[0] = X[0] + c_dt*X[3];
    X_new[1] = X[1] + c_dt*X[4];
    X_new[2] = X[2] + c_dt*X[5];
    X_new[3] = X[3] + c_dt*a[0];
    X_new[4] = X[4] + c_dt*a[1];
    X_new[5] = X[5] + c_dt*a[2];
}

__device__ double calculate_vector_square(double *v){
    return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

__device__ double calculate_dot_product(double *a, double*b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

__device__ void calculate_yoshida_v(double d, double *X_new, double *X, double *a,
                                    double *B,
                                    double eperm, double timestep){
    double d_dt = d*timestep;
    double B_square = calculate_vector_square(B);
    double B_dot_v = calculate_dot_product(B, &X_new[3]);
    double a_new[3];
    get_acceleration_from_lorentz_force(a_new, X_new, B, eperm);
    X_new[3] = X[3] + d_dt*eperm*a_new[0];
    X_new[4] = X[4] + d_dt*eperm*a_new[1];
    X_new[5] = X[5] + d_dt*eperm*a_new[2];
    a_new[0] = a[0] + d_dt*eperm*eperm*(B_dot_v*B[0]-B_square*X[3];
    a_new[1] = a[1] + d_dt*eperm*eperm*(B_dot_v*B[1]-B_square*X[4];
    a_new[2] = a[2] + d_dt*eperm*eperm*(B_dot_v*B[2]-B_square*X[5];
    memcpy(a, a_new, 3*sizeof(double));
}

__device__ void calculate_yoshida_v(double d, double *X_new, double *X, double *a,
                                    double *B, double *E,
                                    double eperm, double timestep){
    double d_dt = d*timestep;
    double B_square = calculate_vector_square(B);
    double B_dot_v = calculate_dot_product(B, &X_new[3]);
    double a_new[3];
    get_acceleration_from_lorentz_force(a_new, X_new, B, E, eperm);
    X_new[3] = X[3] + d_dt*eperm*a_new[0];
    X_new[4] = X[4] + d_dt*eperm*a_new[1];
    X_new[5] = X[5] + d_dt*eperm*a_new[2];
    a_new[0] = a[0] + d_dt*eperm*eperm*(B_dot_v*B[0]-B_square*X[3];
    a_new[1] = a[1] + d_dt*eperm*eperm*(B_dot_v*B[1]-B_square*X[4];
    a_new[2] = a[2] + d_dt*eperm*eperm*(B_dot_v*B[2]-B_square*X[5];
    memcpy(a, a_new, 3*sizeof(double));
}


__device__ void solve_diffeq_by_yoshida(double *X, double *a, double *B,
                                        double eperm, double timestep){
    double X_new[3];
    calculate_yoshida_x(c1, X_new, X, a, timestep);
    calculate_yoshida_v(d1, X_new, X, a, B, eperm, timestep);
    calculate_yoshida_x(c2, X, X_new, a, timestep);
    calculate_yoshida_v(d2, X, X_new, a, B, eperm, timestep);
    calculate_yoshida_x(c3, X_new, X, a, timestep);
    calculate_yoshida_v(d3, X_new, X, a, B, eperm, timestep);
    calculate_yoshida_x(c4, X, X_new, a, timestep);
    memcpy(&X[3], &X_new[3], 3*sizeof(double));
}

__device__ void solve_diffeq_with_efield_by_yoshida(double *X, double *a, double *B, double *E,
                                                    double eperm, double timestep){
    double X_new[3];
    calculate_yoshida_x(c1, X_new, X, a, timestep);
    calculate_yoshida_v(d1, X_new, X, a, B, E, eperm, timestep);
    calculate_yoshida_x(c2, X, X_new, a, timestep);
    calculate_yoshida_v(d2, X, X_new, a, B, E, eperm, timestep);
    calculate_yoshida_x(c3, X_new, X, a, timestep);
    calculate_yoshida_v(d3, X_new, X, a, B, E, eperm, timestep);
    calculate_yoshida_x(c4, X, X_new, a, timestep);
    memcpy(&X[3], &X_new[3], 3*sizeof(double));
}
