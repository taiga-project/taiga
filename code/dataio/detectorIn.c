#include <stdlib.h>
#include <math.h>
#include "../running/detector_postproc.cu"


// set beam inline parameters
void detector_module(double **x_ptr, double *detector, int *detcellid, char *detector_name, int max_blocks, int shot_block_size){
	double *detector_geometry_x;
	double *detector_geometry_y;

	int N_detector_geometry_x = vectorReader(&detector_geometry_x, "input/detector", detector_name, "detx", false);
	int N_detector_geometry_y = vectorReader(&detector_geometry_y, "input/detector", detector_name, "dety", false);

	if (( N_detector_geometry_x > 0) && ( N_detector_geometry_y > 0)) {
		printf("Detector postprocessor module: ON (%d x %d)", N_detector_geometry_x/2, N_detector_geometry_y/2);
		detector_postproc <<< max_blocks, shot_block_size  >>> (x_ptr, detector_geometry_x, N_detector_geometry_x, detector_geometry_y, N_detector_geometry_y, detector, detcellid);
	}else{
		printf("Detector postprocessor module: OFF");
	}
        

}

