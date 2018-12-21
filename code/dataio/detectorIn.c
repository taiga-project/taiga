#include <stdlib.h>
#include <math.h>
#include "../running/detector_postproc.cu"
#include "../running/detector_sum.cu"


// set beam inline parameters
void detector_module(double **x_ptr, double *detector, int *detcellid, char *detector_name, int max_blocks, int shot_block_size, int number_of_particles){
	double *DGX, *dgx;
	double *DGY, *dgy;
	double *DG, *dg;

	int N_dgx = vectorReader(&DGX, "input/detector", detector_name, "detx", false);
	int N_dgy = vectorReader(&DGY, "input/detector", detector_name, "dety", false);
    int N_dg = N_dgx * N_dgy;

	if (( N_dgx > 0) && ( N_dgy > 0)) {
		printf("Detector postprocessor module: ON (%d x %d)", N_dgx/2, N_dgy/2);        
		size_t dimDGX = N_dgx * sizeof(double);     
		size_t dimDGY = N_dgy * sizeof(double);   
		size_t dimDG = N_dg * sizeof(double);
		cudaMalloc((void **) &dgx,  dimDGX);
		cudaMalloc((void **) &dgy,  dimDGY);
		cudaMalloc((void **) &dg,  dimDG);
		cudaMemcpy(dgx, DGX, dimDGX, cudaMemcpyHostToDevice);
		cudaMemcpy(dgy, DGY, dimDGY, cudaMemcpyHostToDevice);        
		detector_postproc <<< max_blocks, shot_block_size  >>> (x_ptr, dgx, N_dgx, dgy, N_dgy, detector, detcellid);
		detector_sum <<<1,1>>>(dg, detcellid, number_of_particles, N_dg);
	}else{
		printf("Detector postprocessor module: OFF");
	}
        

}

