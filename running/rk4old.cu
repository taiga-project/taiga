// Runge--Kutta method

/*__device__ void rk4(double *l_v, double *l_vo, double *l_b, double eperm, double dt){
	
	double a1,a2,a3,a4;
	for(int i=0;i<3;++i){
		a1 = eperm * (  l_vo[(i+1)%3]          * l_b[(i+2)%3] -  l_vo[(i+2)%3]          * l_b[(i+1)%3]);
		a2 = eperm * ( (l_vo[(i+1)%3]+dt/2*a1) * l_b[(i+2)%3] - (l_vo[(i+2)%3]+dt/2*a1) * l_b[(i+1)%3]);
		a3 = eperm * ( (l_vo[(i+1)%3]+dt/2*a2) * l_b[(i+2)%3] - (l_vo[(i+2)%3]+dt/2*a2) * l_b[(i+1)%3]);
		a4 = eperm * ( (l_vo[(i+1)%3]+dt  *a3) * l_b[(i+2)%3] - (l_vo[(i+2)%3]+dt  *a3) * l_b[(i+1)%3]);
		
        l_v[i] = (i+1)%3;//l_vo[i];//+(dt/6)*(a1+2*(a2+a3)+a4);
      

	}
}
*/


// r [0], z [1], t [2]
__device__ double rk4r(double l_vor, double l_voz, double l_vot, double l_br, double l_bz, double l_bt, double eperm){
	
	double l_vr;
	double a1,a2,a3,a4;
	
	a1 = eperm * (  l_voz          * l_bt -  l_vot          * l_bz);
	a2 = eperm * ( (l_voz+dt/2*a1) * l_bt - (l_vot+dt/2*a1) * l_bz);
	a3 = eperm * ( (l_voz+dt/2*a2) * l_bt - (l_vot+dt/2*a2) * l_bz);
	a4 = eperm * ( (l_voz+dt  *a3) * l_bt - (l_vot+dt  *a3) * l_bz);
	
	l_vr = l_vor + (dt/6)*(a1+2*(a2+a3)+a4);
	
	return l_vr;

}

__device__ double rk4z(double l_vor, double l_voz, double l_vot, double l_br, double l_bz, double l_bt, double eperm{
	
	double l_vz;
	double a1,a2,a3,a4;
	
	a1 = eperm * (  l_vot          * l_br -  l_vor          * l_bt);
	a2 = eperm * ( (l_vot+dt/2*a1) * l_br - (l_vor+dt/2*a1) * l_bt);
	a3 = eperm * ( (l_vot+dt/2*a2) * l_br - (l_vor+dt/2*a2) * l_bt);
	a4 = eperm * ( (l_vot+dt  *a3) * l_br - (l_vor+dt  *a3) * l_bt);
	
	l_vz = l_voz + (dt/6)*(a1+2*(a2+a3)+a4);
	
	return l_vz;

}

__device__ double rk4t(double l_vor, double l_voz, double l_vot, double l_br, double l_bz, double l_bt, double eperm){
	
	double l_vt;
	double a1,a2,a3,a4;
	
	a1 = eperm * (  l_vor          * l_bz -  l_voz          * l_br);
	a2 = eperm * ( (l_vor+dt/2*a1) * l_bz - (l_voz+dt/2*a1) * l_br);
	a3 = eperm * ( (l_vor+dt/2*a2) * l_bz - (l_voz+dt/2*a2) * l_br);
	a4 = eperm * ( (l_vor+dt  *a3) * l_bz - (l_voz+dt  *a3) * l_br);
	
	l_vt = l_vot + (dt/6)*(a1+2*(a2+a3)+a4);

	return l_vt;
}


