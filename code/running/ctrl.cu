// Trajectory simulation

//__constant__ int finCounter = 0 ;
__global__ void ctrl(int NR, int NZ, double eperm, double **br_ptr, double **bz_ptr, double **bt_ptr, double **g_ptr, double **x_ptr, double **v_ptr, double *det, int *detcellid, int N_step, double *service_var){
	// thread index
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	detcellid[idx] = -1;
	
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
	
	detcellid[idx] = traj(rg,NR,zg,NZ,sRZT,svRZT,br_ptr,bz_ptr,bt_ptr,eperm,det,N_step,detcellid[idx]);
	//if(idx<20) service_var[idx]=temp;


	x_ptr[0][idx]=sRZT[0];
	x_ptr[1][idx]=sRZT[1];
	x_ptr[2][idx]=sRZT[2];

	v_ptr[0][idx]=svRZT[0];
	v_ptr[1][idx]=svRZT[1];
	v_ptr[2][idx]=svRZT[2];

	if(idx==0){
		 service_var[0] = 42.24;
	}
}

__global__ void ctrl(int NR, int NZ, double eperm, double **br_ptr, double **bz_ptr, double **bt_ptr, double **er_ptr, double **ez_ptr, double **et_ptr, double **g_ptr, double **x_ptr, double **v_ptr, double *det, int *detcellid, int N_step, double *service_var){
	// thread index
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	detcellid[idx] = -1;
	
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
	
	detcellid[idx] = traj(rg,NR,zg,NZ,sRZT,svRZT,br_ptr,bz_ptr,bt_ptr,er_ptr,ez_ptr,et_ptr,eperm,det,N_step,detcellid[idx]);
	//if(idx<20) service_var[idx]=temp;


	x_ptr[0][idx]=sRZT[0];
	x_ptr[1][idx]=sRZT[1];
	x_ptr[2][idx]=sRZT[2];

	v_ptr[0][idx]=svRZT[0];
	v_ptr[1][idx]=svRZT[1];
	v_ptr[2][idx]=svRZT[2];

	if(idx==0){
		 service_var[0] = 42.24;
	}
	
	if(idx==1){
		service_var[1]=er_ptr[0][0];  
	}
	if(idx==2){
		service_var[2]=er_ptr[15][0];  
	}
	if(idx==3){
		service_var[3]=br_ptr[0][0];  
	}
	if(idx==4){
		service_var[4]=br_ptr[15][0];  
	}
}
