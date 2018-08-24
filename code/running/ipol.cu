__device__ int ipol(double l_r, double l_z, double l_t, double l_or, double l_oz, double l_ot, double *det, double *l_x, double l_vr){

	int finished = 0;
	
	double det_R = det[0];
	double det_z = det[1];
	double det_tan = det[3];
	
	double det_dist, det_odist;
	det_dist = (l_r-det_R)+det_tan*(l_z-det_z);
	det_odist = (l_or-det_R)+det_tan*(l_oz-det_z);
	
	if((det_dist*det_odist<=0) && (l_vr>0)){

		double di_rate = (-det_odist)/(det_dist-det_odist);
		l_x[0] = l_or + di_rate*(l_r-l_or);;
		l_x[1] = l_oz + di_rate*(l_z-l_oz);
		l_x[2] = l_ot + di_rate*(l_t-l_ot);

		finished = 1;
	
	}
	
	return finished;
}
