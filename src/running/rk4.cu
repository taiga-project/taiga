// Runge--Kutta method

__device__ void calculate_rk4_coeff(double *X,
                                    double *S, double *S_prev, double rk_weight,
                                    double *B,
                                    double eperm, double timestep){
    for (int i=3; i<6; ++i){
        S[i] = X[i] + rk_weight * S_prev[i];
    }
    
    S[0] = S[3]; // f(R) = vR
    S[1] = S[4]; // f(Z) = vZ
    S[2] = S[5]; // f(T) = vT
    S[3] = eperm*(S[1]*B[2] - S[2]*B[1]); // f(vR) = q/m * (vZ*BT - vT*BZ)
    S[4] = eperm*(S[2]*B[0] - S[0]*B[2]); // f(vZ) = q/m * (vT*BR - vR*BT)
    S[5] = eperm*(S[0]*B[1] - S[1]*B[0]); // f(vT) = q/m * (vR*BZ - vZ*BR)
    
    for (int i=0; i<6; ++i){
        S[i] *= timestep;
    }
}

__device__ void calculate_rk4_coeff(double *X,
                                    double *S, double *S_prev, double rk_weight,
                                    double *B, double *E,
                                    double eperm, double timestep){
    for (int i=3; i<6; ++i){
        S[i] = X[i] + rk_weight * S_prev[i];
    }
    
    S[0] = S[3]; // f(R) = vR
    S[1] = S[4]; // f(Z) = vZ
    S[2] = S[5]; // f(T) = vT
    S[3] = eperm*(E[0] + S[1]*B[2] - S[2]*B[1]); // f(vR) = q/m * (ER + vZ*BT - vT*BZ)
    S[4] = eperm*(E[1] + S[2]*B[0] - S[0]*B[2]); // f(vZ) = q/m * (EZ + vT*BR - vR*BT)
    S[5] = eperm*(E[2] + S[0]*B[1] - S[1]*B[0]); // f(vT) = q/m * (ET + vR*BZ - vZ*BR)
    
    for (int i=0; i<6; ++i){
        S[i] *= timestep;
    }
}

__device__ void solve_diffeq_by_rk4(double *X, double *B,
                                    double eperm, double timestep){
    double S1[6], S2[6], S3[6], S4[6];
    
    calculate_rk4_coeff(X, S1, X,  0.0, B, eperm ,timestep);
    calculate_rk4_coeff(X, S2, S1, 0.5, B, eperm ,timestep);
    calculate_rk4_coeff(X, S3, S2, 0.5, B, eperm ,timestep);
    calculate_rk4_coeff(X, S4, S3, 1.0, B, eperm ,timestep);
    
    for(int i=0; i<6; ++i){
        X[i] = X[i] + (S1[i] + 2*S2[i] + 2*S3[i] + S4[i])/6;
    }
}

__device__ void solve_diffeq_with_efield_by_rk4(double *X, double *B, double *E,
                                                double eperm, double timestep){
    double S1[6], S2[6], S3[6], S4[6];
    
    calculate_rk4_coeff(X, S1, X,  0.0, B, E, eperm ,timestep);
    calculate_rk4_coeff(X, S2, S1, 0.5, B, E, eperm ,timestep);
    calculate_rk4_coeff(X, S3, S2, 0.5, B, E, eperm ,timestep);
    calculate_rk4_coeff(X, S4, S3, 1.0, B, E, eperm ,timestep);
    
    for(int i=0; i<6; ++i){
        X[i] = X[i] + (S1[i] + 2*S2[i] + 2*S3[i] + S4[i])/6;
    }
}
