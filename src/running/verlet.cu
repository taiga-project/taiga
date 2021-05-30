// velocity Verlet
__device__ void solve_diffeq_with_efield_by_verlet(double *X, double *a, double *B,
                                                   double eperm, double timestep){
    double R = X[0] + X[3]*timestep + 0.5*a[0]*timestep*timestep;
    double Z = X[1] + X[4]*timestep + 0.5*a[1]*timestep*timestep;
    double T = X[2] + X[5]*timestep + 0.5*a[2]*timestep*timestep;
    
    double new_a0 = eperm*(X[1]*B[2] - X[2]*B[1]); // f(vR) = q/m * (ER + vZ*BT - vT*BZ)
    double new_a1 = eperm*(X[2]*B[0] - X[0]*B[2]); // f(vZ) = q/m * (EZ + vT*BR - vR*BT)
    double new_a2 = eperm*(X[0]*B[1] - X[1]*B[0]); // f(vT) = q/m * (ET + vR*BZ - vZ*BR)

    X[0] = R;
    X[1] = Z;
    X[2] = T;
    X[3] += (a[0] + new_a0)*timestep*0.5;
    X[4] += (a[1] + new_a1)*timestep*0.5;
    X[5] += (a[2] + new_a2)*timestep*0.5;
    
    a[0] = new_a0;
    a[1] = new_a1;
    a[2] = new_a2;
}

__device__ void solve_diffeq_with_efield_by_verlet(double *X, double *a, double *B, double *E,
                                                   double eperm, double timestep){
    double R = X[0] + X[3]*timestep + 0.5*a[0]*timestep*timestep;
    double Z = X[1] + X[4]*timestep + 0.5*a[1]*timestep*timestep;
    double T = X[2] + X[5]*timestep + 0.5*a[2]*timestep*timestep;
    
    double new_a0 = eperm*(E[0] + X[1]*B[2] - X[2]*B[1]); // f(vR) = q/m * (ER + vZ*BT - vT*BZ)
    double new_a1 = eperm*(E[1] + X[2]*B[0] - X[0]*B[2]); // f(vZ) = q/m * (EZ + vT*BR - vR*BT)
    double new_a2 = eperm*(E[2] + X[0]*B[1] - X[1]*B[0]); // f(vT) = q/m * (ET + vR*BZ - vZ*BR)

    X[0] = R;
    X[1] = Z;
    X[2] = T;
    X[3] += (a[0] + new_a0)*timestep*0.5;
    X[4] += (a[1] + new_a1)*timestep*0.5;
    X[5] += (a[2] + new_a2)*timestep*0.5;
    
    a[0] = new_a0;
    a[1] = new_a1;
    a[2] = new_a2;
}
