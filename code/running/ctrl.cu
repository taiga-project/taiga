// Trajectory simulation

//__constant__ int finCounter = 0 ;

__global__ void banCtrl(int NR, int NZ, double **br_ptr, double **bz_ptr, double **bt_ptr, double **g_ptr, double **x_ptr,double *bd1, double *bd2){

	// thread index
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//		double temp;
	
	// grid pointer
	double *rg, *zg;
	rg = g_ptr[0];	
	zg = g_ptr[1];
	
	//double valR=0.72, valZ=0.0;
	double sRZT[3], l_b[3];
	//double *trp;
	//int i,j;

	sRZT[0] = x_ptr[0][idx];
	sRZT[1] = x_ptr[1][idx];
	sRZT[2] = x_ptr[2][idx];

	
	banTraj(rg,NR,zg,NZ,sRZT,br_ptr,bz_ptr,bt_ptr,l_b);

	bd1[idx]=l_b[1]/l_b[0];
	bd2[idx]=l_b[2]/l_b[0];



}

__global__ void ctrl(int NR, int NZ, double **br_ptr, double **bz_ptr, double **bt_ptr, double **g_ptr, double **x_ptr, double **v_ptr, double *tmp, double eperm, double l_ri, int N_step){
	// thread index
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//		double temp;
	
	// grid pointer
	double *rg, *zg;
	rg = g_ptr[0];	
	zg = g_ptr[1];
	
	//double valR=0.72, valZ=0.0;
	double sRZT[3], svRZT[3];
	//double *trp;
	//int i,j;

	sRZT[0] = x_ptr[0][idx];
	sRZT[1] = x_ptr[1][idx];
	sRZT[2] = x_ptr[2][idx];

	svRZT[0] = v_ptr[0][idx];
	svRZT[1] = v_ptr[1][idx];
	svRZT[2] = v_ptr[2][idx];
	
	traj(rg,NR,zg,NZ,sRZT,svRZT,br_ptr,bz_ptr,bt_ptr,eperm,l_ri,N_step);
	//if(idx<20) tmp[idx]=temp;


	x_ptr[0][idx]=sRZT[0];
	x_ptr[1][idx]=sRZT[1];
	x_ptr[2][idx]=sRZT[2];

	v_ptr[0][idx]=svRZT[0];
	v_ptr[1][idx]=svRZT[1];
	v_ptr[2][idx]=svRZT[2];

	if(idx==0){
		 tmp[0] = 42.24;
	}
}

__global__ void ctrl(int NR, int NZ, double **br_ptr, double **bz_ptr, double **bt_ptr, double **er_ptr, double **ez_ptr, double **et_ptr, double **g_ptr, double **x_ptr, double **v_ptr, double *tmp, double eperm, double l_ri, int N_step){
	// thread index
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//		double temp;
	
	// grid pointer
	double *rg, *zg;
	rg = g_ptr[0];	
	zg = g_ptr[1];
	
	//double valR=0.72, valZ=0.0;
	double sRZT[3], svRZT[3];
	//double *trp;
	//int i,j;

	sRZT[0] = x_ptr[0][idx];
	sRZT[1] = x_ptr[1][idx];
	sRZT[2] = x_ptr[2][idx];

	svRZT[0] = v_ptr[0][idx];
	svRZT[1] = v_ptr[1][idx];
	svRZT[2] = v_ptr[2][idx];
	
	traj(rg,NR,zg,NZ,sRZT,svRZT,br_ptr,bz_ptr,bt_ptr,er_ptr,ez_ptr,et_ptr,eperm,l_ri,N_step);
	//if(idx<20) tmp[idx]=temp;


	x_ptr[0][idx]=sRZT[0];
	x_ptr[1][idx]=sRZT[1];
	x_ptr[2][idx]=sRZT[2];

	v_ptr[0][idx]=svRZT[0];
	v_ptr[1][idx]=svRZT[1];
	v_ptr[2][idx]=svRZT[2];

	if(idx==0){
		 tmp[0] = 42.24;
	}
}

/*
	sRZT[0]=valR;
	sRZT[1]=valZ;
	sRZT[2]=0;
	
//	double gyok2 = 0.7071;
	
	//br2[idx]=br1[idx]*2;
	//b_ptr[1][idx]=b_ptr[0][idx]/5;

	if(idx==0) tmp[0] = 42.424242;
	if(idx==1){
		 tmp[1] = traj(rg,NR,zg,NZ,sRZT,svRZT,br_ptr,bz_ptr,bt_ptr,eperm);
		 tmp[20]=sRZT[0];
 		 tmp[21]=sRZT[1];
 		 tmp[22]=sRZT[2];
	 }
	//if(idx==2) tmp[2] = 2014.1008;
	if(idx==2) tmp[2] = g_ptr[0][0];
	
	
	if(idx==3){

		for(i=0;(rg[i]<valR)&&(i<NR-1);i++){;}
		
		tmp[3] = rg[i-1];
		tmp[5] = rg[i];//(double)i;
		tmp[6] = (double)i;
		
	}
	
	if(idx==4) tmp[4] = valR;
	if(idx==7){

		for(j=0;(zg[j]<valZ)&&(j<NZ-1);j++){;}
		tmp[7] = zg[j-1];
		tmp[9] = zg[j];//(double)i;
		tmp[10] = (double)j;
	}
	if(idx==8) tmp[8] = valZ;
	
	if(idx==11) tmp[11] = svRZT[0];
	if(idx==12) tmp[12] = svRZT[1];
	if(idx==13) tmp[13] = svRZT[2];
	
	
	
	if(idx==14) tmp[14] = (double)(NZ-1);

	if(idx==18) tmp[18] = NR;
	
	if(idx==19) tmp[19] = NZ;
	
//	if(idx==0) tmp[0] = idx+0.42;*/
//	return 1;	
