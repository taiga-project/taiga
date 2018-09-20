#include <stdlib.h>
#include <math.h>
#include "../running/detector_postproc.cu"


// set beam inline parameters
void detector_module(double **x_ptr, double *detector, int *detcellid, char *detector_name, int max_blocks, int shot_block_size){
	double *DGX, *dgx;
	double *DGY, *dgy;

	int N_dgx = vectorReader(&DGX, "input/detector", detector_name, "detx", false);
	int N_dgy = vectorReader(&DGY, "input/detector", detector_name, "dety", false);    

	if (( N_dgx > 0) && ( N_dgy > 0)) {
		printf("Detector postprocessor module: ON (%d x %d)", N_dgx/2, N_dgy/2);        
		size_t dimDGX = N_dgx * sizeof(double);     
		size_t dimDGY = N_dgy * sizeof(double);
		cudaMalloc((void **) &dgx,  dimDGX);
		cudaMalloc((void **) &dgy,  dimDGY);
		cudaMemcpy(dgx, DGX, dimDGX, cudaMemcpyHostToDevice);
		cudaMemcpy(dgy, DGY, dimDGY, cudaMemcpyHostToDevice);        
		detector_postproc <<< max_blocks, shot_block_size  >>> (x_ptr, dgx, N_dgx, dgy, N_dgy, detector, detcellid);
	}else{
		printf("Detector postprocessor module: OFF");
	}
        

}

