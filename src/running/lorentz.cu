__device__ void get_acceleration_from_lorentz_force(double *a, double *v,
                                              double *B, double *E,
                                              double eperm){
    a[0] = eperm*(E[0] + v[1]*B[2] - v[2]*B[1]); // f(vR) = q/m * (ER + vZ*BT - vT*BZ)
    a[1] = eperm*(E[1] + v[2]*B[0] - v[0]*B[2]); // f(vZ) = q/m * (EZ + vT*BR - vR*BT)
    a[2] = eperm*(E[2] + v[0]*B[1] - v[1]*B[0]); // f(vT) = q/m * (ET + vR*BZ - vZ*BR)
}

__device__ void get_acceleration_from_lorentz_force(double *a, double *v,
                                              double *B,
                                              double eperm){
    a[0] = eperm*(v[1]*B[2] - v[2]*B[1]); // f(vR) = q/m * (vZ*BT - vT*BZ)
    a[1] = eperm*(v[2]*B[0] - v[0]*B[2]); // f(vZ) = q/m * (vT*BR - vR*BT)
    a[2] = eperm*(v[0]*B[1] - v[1]*B[0]); // f(vT) = q/m * (vR*BZ - vZ*BR)
}

