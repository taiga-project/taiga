
//double gyok2;

__device__ /*double int*/ void copyLocal(double *rg, int NR, double *zg, int NZ, double l_r, double l_z, int *rzci, double *lp_br, double *lp_bz, double *lp_bt,  double **br_ptr, double **bz_ptr, double **bt_ptr){
	int rci, zci;
	int i, i2;
	
	for(rci=0;(rg[rci+1]<l_r)&&(rci<NR-1);rci++){;}	
	
	for(zci=0;(zg[zci+1]<l_z)&&(zci<NR-1);zci++){;}	
	
	// Particle leave out the cell
	if ((rzci[0]!=rci)||(rzci[1]!=zci)){
		rzci[0]=rci; 	
		rzci[1]=zci; 	
	
		for(i=0;i<16;i++){
	
			i2 = (rzci[0])*(NZ-1)+rzci[1];
		
			lp_br[i]=br_ptr[i][i2];			
			lp_bz[i]=bz_ptr[i][i2];			
			lp_bt[i]=bt_ptr[i][i2];
		
		}
	}
	
	//return i2;	
}

__device__ /*double int*/ void copyLocal(double *rg, int NR, double *zg, int NZ, double l_r, double l_z, int *rzci, double *lp_br, double *lp_bz, double *lp_bt,  double **br_ptr, double **bz_ptr, double **bt_ptr,
										double *lp_er, double *lp_ez, double *lp_et,  double **er_ptr, double **ez_ptr, double **et_ptr){
	int rci, zci;
	int i, i2;
	
	for(rci=0;(rg[rci+1]<l_r)&&(rci<NR-1);rci++){;}	
	
	for(zci=0;(zg[zci+1]<l_z)&&(zci<NR-1);zci++){;}	
	
	// Particle leave out the cell
	if ((rzci[0]!=rci)||(rzci[1]!=zci)){
		rzci[0]=rci; 	
		rzci[1]=zci; 	
	
		for(i=0;i<16;i++){
	
			i2 = (rzci[0])*(NZ-1)+rzci[1];
		
			lp_br[i]=br_ptr[i][i2];			
			lp_bz[i]=bz_ptr[i][i2];			
			lp_bt[i]=bt_ptr[i][i2];
			lp_er[i]=er_ptr[i][i2];			
			lp_ez[i]=ez_ptr[i][i2];			
			lp_et[i]=et_ptr[i][i2];
		
		}
	}
	
	//return i2;	
}


__device__ double localField(double *lp_b, double dr, double dz){

/* MATLAB CODE:
    sample2(3) =c11(bs1,bs2)*dsx^3*dsy^3 + c12(bs1,bs2)*dsx^3*dsy^2 + c13(bs1,bs2)*dsx^3*dsy + c14(bs1,bs2)*dsx^3 + ...
                c21(bs1,bs2)*dsx^2*dsy^3 + c22(bs1,bs2)*dsx^3*dsy^2 + c23(bs1,bs2)*dsx^2*dsy + c24(bs1,bs2)*dsx^2 + ...
                c31(bs1,bs2)*dsx  *dsy^3 + c32(bs1,bs2)*dsx  *dsy^2 + c33(bs1,bs2)*dsx  *dsy + c34(bs1,bs2)*dsx    + ...
                c41(bs1,bs2)      *dsy^3 + c42(bs1,bs2)      *dsy^2 + c43(bs1,bs2)      *dsy + c44(bs1,bs2);*/
                
                
    double blocal = 0.0, tmp[16] ;
    for(int i=0;i<4;i++){
 	    for(int j=0;j<4;j++){
		    tmp[i*4+j] = lp_b[i*4+j]*pow(dr,3-i)*pow(dz,3-j);
			
        }
    }   
    
   	for(int i=0;i<4;i++){
    	for(int j=0;j<4;j++){
		    blocal+=tmp[i*4+j];
        }
    }   
    
	return blocal;
}

/*

l_x local coordinates

*/
__device__ void traj(double *rg, int NR, double *zg, int NZ, double *l_x, double *l_v, double **br_ptr, double **bz_ptr, double **bt_ptr, double eperm, double *det, int N_step){

	// next grid
	int rzci[2];
	rzci[0]=-1;
	rzci[1]=-1;
		
	double lp_br[16];
	double lp_bz[16];
	double lp_bt[16];	
	
	double l_br=0,l_bz,l_bt;
	double dr,dz;
//	double drT;
	double l_rT;
	
	double l_vr, l_vz, l_vt, l_vor, l_voz, l_vot;
	double l_r,  l_z,   l_t, l_or,  l_oz,  l_ot;

	// X = [R,vR,z,vZ,T,vT]
	double X[6];
	
	int finished = 0;
	
	l_r  = l_x[0];
	l_z  = l_x[1];
	l_t  = l_x[2];
	
	l_vr = l_v[0];
	l_vz = l_v[1];
	l_vt = l_v[2];

	l_x[0]=0.0;
	
	int loopi;
	for (loopi=0;(loopi<N_step && (!finished));loopi++){
		// Get local magnetic field	

		l_rT = cyl2tor_coord(l_r, l_t);
		copyLocal(rg,NR,zg,NZ,l_rT,l_z,rzci,lp_br,lp_bz,lp_bt,br_ptr,bz_ptr,bt_ptr);	
		
		//dr = l_r-rg[rzci[0]];
		dr = l_rT-rg[rzci[0]];
		dz = l_z-zg[rzci[1]];
	
		l_br =  localField(lp_br,dr,dz);
		l_bz =  localField(lp_bz,dr,dz);
		l_bt =  localField(lp_bt,dr,dz);

		l_br = cyl2tor_rad(l_br, l_bt, l_r, l_t);
		l_bt = cyl2tor_tor(l_br, l_bt, l_r, l_t);
	

		// archivate coordinates
		l_or  = l_r;	l_oz  = l_z;	l_ot  = l_t;
		l_vor = l_vr;	l_voz = l_vz;	l_vot = l_vt;

		// new coordinates
		X[0] = l_r;
		X[2] = l_z;
		X[4] = l_t;

		// new speed components
		X[1] = l_vr;
		X[3] = l_vz;
		X[5] = l_vt;
	
		rk4lin(X, l_br, l_bz, l_bt, eperm);
		
		// new coordinates
		l_r = X[0];
		l_z = X[2];
		l_t = X[4];

		// new speed components
		l_vr = X[1];
		l_vz = X[3];
		l_vt = X[5];
	
		// finished? (interpolation)
		finished = ipol(l_r, l_z, l_t, l_or, l_oz, l_ot, det, l_x, l_vr);

	}
	
	l_v[0] = l_vr;
	l_v[1] = l_vz;
	l_v[2] = l_vt;

	if (!finished){
		l_x[0] = NAN;//l_r;
		l_x[1] = NAN;//l_z;
		l_x[2] = NAN;//l_t;
	}
	
}

__device__ void traj(double *rg, int NR, double *zg, int NZ, double *l_x, double *l_v, double **br_ptr, double **bz_ptr, double **bt_ptr, double **er_ptr, double **ez_ptr, double **et_ptr, double eperm, double *det, int N_step){

	// next grid
	int rzci[2];
	rzci[0]=-1;
	rzci[1]=-1;
		
	double lp_br[16];
	double lp_bz[16];
	double lp_bt[16];

	double lp_er[16];
	double lp_ez[16];
	double lp_et[16];	
	
	double l_br=0,l_bz=0,l_bt=0;
	double l_er=0,l_ez=0,l_et=0;
	double dr,dz;
	double l_rT;
	
	double l_vr, l_vz, l_vt, l_vor, l_voz, l_vot;
	double l_r,  l_z,   l_t, l_or,  l_oz,  l_ot;

	double X[6];
	
	int finished = 0;		

	l_r  = l_x[0];
	l_z  = l_x[1];
	l_t  = l_x[2];
	
	l_vr = l_v[0];
	l_vz = l_v[1];
	l_vt = l_v[2];

	l_x[0]=0.0;
	
	int loopi;
	for (loopi=0;(loopi<N_step && (!finished));loopi++){
		// Get local magnetic field	

		l_rT = cyl2tor_coord(l_r, l_t);
		copyLocal(rg,NR,zg,NZ,l_rT,l_z,rzci,lp_br,lp_bz,lp_bt,br_ptr,bz_ptr,bt_ptr,lp_er,lp_ez,lp_et,er_ptr,ez_ptr,et_ptr);	
		
		dr = l_rT-rg[rzci[0]];
		dz = l_z-zg[rzci[1]];
	
		l_br =  localField(lp_br,dr,dz);
		l_bz =  localField(lp_bz,dr,dz);
		l_bt =  localField(lp_bt,dr,dz);
		l_br = cyl2tor_rad(l_br, l_bt, l_r, l_t);
		l_bt = cyl2tor_tor(l_br, l_bt, l_r, l_t);
		
		l_er =  localField(lp_er,dr,dz);
		l_ez =  localField(lp_ez,dr,dz);
		l_et =  localField(lp_et,dr,dz);		
		l_er = cyl2tor_rad(l_er, l_et, l_r, l_t);
		l_et = cyl2tor_tor(l_er, l_et, l_r, l_t);

		// archivate coordinates
		l_or  = l_r;	l_oz  = l_z;	l_ot  = l_t;
		l_vor = l_vr;	l_voz = l_vz;	l_vot = l_vt;

		// new coordinates
		X[0] = l_r;
		X[2] = l_z;
		X[4] = l_t;

		// new speed components
		X[1] = l_vr;
		X[3] = l_vz;
		X[5] = l_vt;
	
		rk4lin(X, l_br, l_bz, l_bt, l_er, l_ez, l_et, eperm);
		
		// new coordinates
		l_r = X[0];
		l_z = X[2];
		l_t = X[4];

		// new speed components
		l_vr = X[1];
		l_vz = X[3];
		l_vt = X[5];
	
	
		// finished? (interpolation)
		finished = ipol(l_r, l_z, l_t, l_or, l_oz, l_ot, det, l_x, l_vr);

	}

	
	l_v[0] = l_vr;
	l_v[1] = l_vz;
	l_v[2] = l_vt;
	

	if (!finished){
		l_x[0] = NAN;//l_r;
		l_x[1] = NAN;//l_z;
		l_x[2] = NAN;//l_t;
	}
	
}

__device__ void banTraj(double *rg, int NR, double *zg, int NZ, double *l_x,  double **br_ptr, double **bz_ptr, double **bt_ptr,  double *l_b){

	// next grid
	int rzci[2];
	rzci[0]=-1;
	rzci[1]=-1;
		
	double lp_br[16];
	double lp_bz[16];
	double lp_bt[16];	
	
	double l_br=0,l_bz,l_bt;
	double dr,dz;
//	double drT;
	double l_rT;
	
	double l_r,  l_z,   l_t;

	// X = [R,vR,z,vZ,T,vT]
	double X[6];
		
	

	l_r  = l_x[0];
	l_z  = l_x[1];
	l_t  = l_x[2];
	

	l_x[0]=0.0;
	

	l_rT = cyl2tor_coord(l_r, l_t);
	copyLocal(rg,NR,zg,NZ,l_rT,l_z,rzci,lp_br,lp_bz,lp_bt,br_ptr,bz_ptr,bt_ptr);	
	
	//dr = l_r-rg[rzci[0]];
	dr = l_rT-rg[rzci[0]];
	dz = l_z-zg[rzci[1]];

	l_br =  localField(lp_br,dr,dz);
	l_bz =  localField(lp_bz,dr,dz);
	l_bt =  localField(lp_bt,dr,dz);

	l_br = cyl2tor_rad(l_br, l_bt, l_r, l_t);
	l_bt = cyl2tor_tor(l_br, l_bt, l_r, l_t);
	
	l_b[0] = l_br;
	l_b[1] = l_bz;
	l_b[2] = l_bt;
		


	
}
