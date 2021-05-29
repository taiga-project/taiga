// Runge--Kutta method

__device__ void calculate_rk4_coeff(double *X,
                                    double *S, double *S_prev, double rk_weight,
                                    double l_br, double l_bz,  double l_bt,
                                    double eperm, double timestep){
    for (int i=3; i<6; ++i){
        S[i] = X[i] + rk_weight * S_prev[i];
    }
    
    S[0] = S[3]; // f(R) = vR
    S[1] = S[4]; // f(Z) = vZ
    S[2] = S[5]; // f(T) = vT
    S[3] = eperm*(S[1]*l_bt - S[2]*l_bz); // f(vR) = q/m * (vZ*BT - vT*BZ)
    S[4] = eperm*(S[2]*l_br - S[0]*l_bt); // f(vZ) = q/m * (vT*BR - vR*BT)
    S[5] = eperm*(S[0]*l_bz - S[1]*l_br); // f(vT) = q/m * (vR*BZ - vZ*BR)
    
    for (int i=0; i<6; ++i){
        S[i] *= timestep;
    }
}

__device__ void calculate_rk4_coeff(double *X,
                                    double *S, double *S_prev, double rk_weight,
                                    double l_br, double l_bz,  double l_bt,
                                    double l_er, double l_ez, double l_et,
                                    double eperm, double timestep){
    for (int i=3; i<6; ++i){
        S[i] = X[i] + rk_weight * S_prev[i];
    }
    
    S[0] = S[3]; // f(R) = vR
    S[1] = S[4]; // f(Z) = vZ
    S[2] = S[5]; // f(T) = vT
    S[3] = eperm*(l_er + S[1]*l_bt - S[2]*l_bz); // f(vR) = q/m * (ER + vZ*BT - vT*BZ)
    S[4] = eperm*(l_ez + S[2]*l_br - S[0]*l_bt); // f(vZ) = q/m * (EZ + vT*BR - vR*BT)
    S[5] = eperm*(l_et + S[0]*l_bz - S[1]*l_br); // f(vT) = q/m * (ET + vR*BZ - vZ*BR)
    
    for (int i=0; i<6; ++i){
        S[i] *= timestep;
    }
}

__device__ void solve_diffeq_by_rk4(double *X,
                                    double l_br, double l_bz, double l_bt,
                                    double eperm, double timestep){
    double S1[6], S2[6], S3[6], S4[6];
    
    calculate_rk4_coeff(X, S1, X,  0.0, l_br, l_bz, l_bt, eperm ,timestep);
    calculate_rk4_coeff(X, S2, S1, 0.5, l_br, l_bz, l_bt, eperm ,timestep);
    calculate_rk4_coeff(X, S3, S2, 0.5, l_br, l_bz, l_bt, eperm ,timestep);
    calculate_rk4_coeff(X, S4, S3, 1.0, l_br, l_bz, l_bt, eperm ,timestep);
    
    for(int i=0; i<6; ++i){
        X[i] = X[i] + (S1[i] + 2*S2[i] + 2*S3[i] + S4[i])/6;
    }
}

__device__ void solve_diffeq_with_efield_by_rk4(double *X,
                                                double l_br, double l_bz, double l_bt,
                                                double l_er, double l_ez, double l_et,
                                                double eperm, double timestep){
    double S1[6], S2[6], S3[6], S4[6];
    
    calculate_rk4_coeff(X, S1, X,  0.0, l_br, l_bz, l_bt, l_er, l_ez, l_et, eperm ,timestep);
    calculate_rk4_coeff(X, S2, S1, 0.5, l_br, l_bz, l_bt, l_er, l_ez, l_et, eperm ,timestep);
    calculate_rk4_coeff(X, S3, S2, 0.5, l_br, l_bz, l_bt, l_er, l_ez, l_et, eperm ,timestep);
    calculate_rk4_coeff(X, S4, S3, 1.0, l_br, l_bz, l_bt, l_er, l_ez, l_et, eperm ,timestep);
    
    for(int i=0; i<6; ++i){
        X[i] = X[i] + (S1[i] + 2*S2[i] + 2*S3[i] + S4[i])/6;
    }
}
