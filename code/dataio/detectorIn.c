#include <stdlib.h>
#include <math.h>


// set beam inline parameters
void detector_module(double **x_ptr, double *detector, int *detcellid, char *detector_name){
	double *detector_geometry_x;
	double *detector_geometry_y;

	int N_detector_geometry_x = vectorReader(&detector_geometry_x, "input/detector", detector_name, "detx", false);
	int N_detector_geometry_y = vectorReader(&detector_geometry_y, "input/detector", detector_name, "dety", false);

	detector_postproc <<< max_blocks, shot.block_size  >>> (x_ptr, detector_geometry_x, N_detector_geometry_x, detector_geometry_y, N_detector_geometry_y, detector, detcellid);
	

}

