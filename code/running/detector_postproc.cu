// Trajectory simulation

//__constant__ int finCounter = 0 ;
__global__ void detector_postproc(double **x_ptr, double *det_x, int N_det_x, double *det_y, int N_det_y, double *det, int *detcellid/*, double *service_var*/){
	// thread index
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (detcellid[idx] == 0){
		double x,y;
		
		x = x_ptr[2][idx] - det[2];                x*=100;
		y = (x_ptr[1][idx] - det[1]) / det[3];     y*=100;    
		
		int x_cellid = -1, y_cellid = -1;
		
		//if ((x >= det_x[0]) && (x <= det_x[N_det_x-1]) && (y >= det_y[0]) && (y <= det_y[N_det_y-1])) {
			for (int i=0; i<N_det_x/2; i++) {
				if ((x >= det_x[2*i]) && (x <= det_x[2*i+1]))	x_cellid = i;
			}
			for (int j=0; j<N_det_y/2; j++) {
				if ((y >= det_y[2*j]) && (y <= det_y[2*j+1]))	y_cellid = j;
			}
		//}
		
		if ((x_cellid >= 0) && (y_cellid >= 0)) {
			detcellid[idx] = x_cellid*1000+y_cellid;//(x_cellid-1)*N_det_y/2 + y_cellid;
		}
	}
	/*if(idx==0){
		 service_var[0] = 42.24;
	}*/
}
//https://www.reverbnation.com/riky9550
