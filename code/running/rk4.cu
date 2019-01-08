// Runge--Kutta method

/*
Linearized operator of RK4

-->
	*X:		previous coordinates (input of rk4lin)
	*S:		output of rk4op (RK4 component)
	*Sp:		previous RK4 coefficient 
	mSp:		multiplier for X argument
	l_br:		radial component of local poloidal field
	l_bz:		vertical component of local poloidal field
	l_bt:		local toroidal field
	eperm:		q/m
	dt:		timestep

<--
	*S:		f(X  + mSp * Sp) * dt

*/
__device__ void rk4op(double *X, double *S, double *Sp, double mSp, double l_br, double l_bz, double l_bt, double eperm){

	for (int i=0;i<6;i++){
		S[i] = X[i] + mSp * Sp[i];
	}
	
	S[0] = S[1]; // f(R) = vR
	S[2] = S[3]; // f(Z) = vZ
	S[4] = S[5]; // f(T) = vT

	S[1] = eperm*(S[2]*l_bt - S[4]*l_bz);		// q/m * (X4*B3-X6*B2)  = q/m * (v2*B3-v3*B2)
	S[3] = eperm*(S[4]*l_br - S[0]*l_bt);		// q/m * (X6*B1-X2*B3)  = q/m * (v3*B1-v1*B3)
	S[5] = eperm*(S[0]*l_bz - S[2]*l_br);		// q/m * (X2*B2-X4*B1)  = q/m * (v1*B2-v2*B1)

	for (int i=0;i<6;i++){
		S[i] *= dt;
	}
}

__device__ void rk4op(double *X, double *S, double *Sp, double mSp, double l_br, double l_bz,  double l_bt, double l_er, double l_ez, double l_et, double eperm){

	for (int i=0;i<6;i++){
		S[i] = X[i] + mSp * Sp[i];
	}
	
	S[0] = S[1]; // f(R) = vR
	S[2] = S[3]; // f(Z) = vZ
	S[4] = S[5]; // f(T) = vT
	S[1] = eperm*(l_er + S[2]*l_bt - S[4]*l_bz);		// q/m * (E1+X4*B3-X6*B2)  = q/m * (E1+v2*B3-v3*B2)
	S[3] = eperm*(l_ez + S[4]*l_br - S[0]*l_bt);		// q/m * (E2+X6*B1-X2*B3)  = q/m * (E2+v3*B1-v1*B3)
	S[5] = eperm*(l_et + S[0]*l_bz - S[2]*l_br);		// q/m * (E3+X2*B2-X4*B1)  = q/m * (E3+v1*B2-v2*B1)

	for (int i=0;i<6;i++){
		S[i] *= dt;
	}
}

__device__ void rk4lin(double *X, double l_br, double l_bz, double l_bt, double eperm){
	
	double S1[6],S2[6],S3[6],S4[6];

	rk4op(X, S1, X,  0.0, l_br, l_bz, l_bt, eperm);
	rk4op(X, S2, S1, 0.5, l_br, l_bz, l_bt, eperm);
	rk4op(X, S3, S2, 0.5, l_br, l_bz, l_bt, eperm);
	rk4op(X, S4, S3, 1.0, l_br, l_bz, l_bt, eperm);

	for(int i=0;i<6;i++){
		X[i] = X[i] + (S1[i] + 2*S2[i] + 2*S3[i] + S4[i])/6;		
	}
}

__device__ void rk4lin(double *X, double l_br, double l_bz, double l_bt, double l_er, double l_ez, double l_et, double eperm){
	
	double S1[6],S2[6],S3[6],S4[6];

	rk4op(X, S1, X,  0.0, l_br, l_bz, l_bt, l_er, l_ez, l_et, eperm);
	rk4op(X, S2, S1, 0.5, l_br, l_bz, l_bt, l_er, l_ez, l_et, eperm);
	rk4op(X, S3, S2, 0.5, l_br, l_bz, l_bt, l_er, l_ez, l_et, eperm);
	rk4op(X, S4, S3, 1.0, l_br, l_bz, l_bt, l_er, l_ez, l_et, eperm);

	for(int i=0;i<6;i++){
		X[i] = X[i] + (S1[i] + 2*S2[i] + 2*S3[i] + S4[i])/6;		
	}
}
