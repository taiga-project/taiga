__device__ double glvar  = 3.14159265358979324;

//#define LSS 441
#define LS 2	
const int LS2 = 2*LS+1;
const int LSS = LS2*LS2;
//double gyok2;

__device__ double copyLocal(double *rg, int NR, double *zg, int NZ, double l_r, double l_z, int *rzci){
	int rci, zci;
	
	for(rci=0;(rg[rci]<l_r)&&(rci<NR-1);rci++){;}
	
	if(rci<LS)    rci = LS;      // left   border
	if(NR-rci<LS) rci = NR-LS-1; // right  border
	
	
	for(zci=0;(zg[zci]<l_z)&&(zci<NR-1);zci++){;}
	
	if(zci<LS)    zci = LS;      // bottom border
	if(NZ-zci<LS) zci = NZ-LS-1; // top    border
	
	rzci[0]=rci-LS;  /* left edge of neighbours */
	rzci[1]=zci-LS;  /* top edge of neighbours */
		
	return rci+zci/100.0;
	//zg[zci];	
}


__device__ double fill(double *rg, int NR, double *zg, int NZ, double l_r, double l_z, double **br_ptr, double **bz_ptr, double **bt_ptr){
	//double //l_br[100];
	double tmp;
	
	double rl[LS2+1]; /*first and last ==6 (LS2==5)*/
	double zl[LS2+1];
	int rzci[2];
	

	
	
	double l_br[16][LSS];
	double l_bz[16][LSS];	
	double l_bt[16][LSS];
	
	
	//return tmp;
	copyLocal(rg,NR,zg,NZ,l_r,l_z,rzci);
	
	
	tmp=0;
	
	/*for(int r_i=rzci[0]-2;r_i<=rzci[0]+2;r_i++){
	
		for(int z_i=rzci[1]-2;z_i<=rzci[1]+2;z_i++){
			if (r_i==rzci[0]){
				tmp*=100;
				tmp+=z_i;
			}
		}
	}*/
	
	
	int r_i,z_i,i,i2,i2l;
	for (r_i=0;r_i<LS2;r_i++){
		for (z_i=0;z_i<LS2;z_i++){
			for(i=0;i<16;i++){
			
				i2 = (r_i+rzci[0])*(NZ-1)+z_i+rzci[1];//i2=8812;
				i2l = r_i*LS2+z_i;
				
				l_br[i][i2l]=br_ptr[i][i2];			
				l_bz[i][i2l]=bz_ptr[i][i2];			
				l_bt[i][i2l]=bt_ptr[i][i2];
				
			}
		}
	}
	
	for (i=0;i<=LS2;i++){
		rl[i] = rg[rzci[0]+i];		
		zl[i] = zg[rzci[1]+i];
	}
	
	
	
	
	
	
	
	
	
	// END OF FILL
	
	// BEGIN OF STEP
	
	
	double dr,dz,dtemp;

	dtemp=0;
	
	for(r_i=0;(r_i<LS2)&&(dtemp>=0);r_i++){
		dr = dtemp;
		dtemp = l_r-rl[r_i]; 	
	}
	
	r_i--;
	
//	r_i=i-1;
	dtemp=0;
	
	for(z_i=0;(z_i<LS2)&&(dtemp>=0);z_i++){
		dz = dtemp;
		dtemp = l_z-zl[z_i]; 	
	}
	z_i--;
	
	
	
	//END OF STEP
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
/*	for(r_i=0;(r_i<5)&&(<l_r);r_i++){
	;
	}
		*/
	
	//return(LS+rzci[0])*(NZ-1)+LS+rzci[1];
	return l_br[0][LS2*LS+LS];
	
	
//	return LS2*LS+LS;
	
	//return 123456789.012345;
//	return z_i;	
//	return LS2;

//	return br_ptr[0][0];	
	
	
	//	return (double)rzci[1];
}
