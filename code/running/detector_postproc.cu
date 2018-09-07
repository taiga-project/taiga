// Trajectory simulation

//__constant__ int finCounter = 0 ;
__global__ void detector_postproc(double **x_ptr, double *det, int *detcellid, double *service_var){
	// thread index
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	// grid pointer
	double *tg, *zg;
	tg = g_ptr[0];	
	zg = g_ptr[1];
	
	double sRZT[3];

	coordinate[0] = x_ptr[0][idx];
	coordinate[1] = x_ptr[1][idx];
	coordinate[2] = x_ptr[2][idx];

	if(idx==0){
		 service_var[0] = 42.24;
	}
}
