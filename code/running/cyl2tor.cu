/* rad, z, tor */


/*__device__ void cyl2tor(double *fieldC, double *fieldT, double xtC, double xrC){


	double tanxr = xtC/xrC;
	double cosxr = 1/sqrt(1+tanxr*tanxr);
	// sinxr = tanxr*cosxr


	// radT = radC * cos(phi)+  torC  * sin(phi)
	fieldT[0]  = fieldC[0] * cosxr + fieldC[2] * tanxr*cosxr;

	// zT   = zC
	fieldT[1]  = fieldC[1];

	// torT = -radC  *   sin(phi)   + torC   * cos(phi)
	fieldT[2]  = -fieldC[0] * tanxr*cosxr + fieldC[2] * cosxr;
	
}*/

// field: B_R{CYL} -> B_R{TOR}
__device__ double cyl2tor_rad(double l_br, double l_bt, double l_r, double l_t){

	double tanxr = l_t/l_r;
//	double cosxr = (1+tanxr*tanxr);
	double cosxr = 1/sqrt(1+tanxr*tanxr);





	// radT = radC * cos(phi) + torC * sin(phi)
//	return    l_br * cosxr    + l_bt * tanxr*cosxr;	// ez a régi
	return    l_br * cosxr    - l_bt * tanxr*cosxr;
}


// field: B_tor{CYL} -> B_tor{TOR}
__device__ double cyl2tor_field(double l_br, double l_bt, double l_r, double l_t){

	double tanxr = l_t/l_r;
//	double cosxr = (1+tanxr*tanxr);
	double cosxr = 1/sqrt(1+tanxr*tanxr);

	// torT = -radC *   sin(phi)  + torC * cos(phi)
//	return    -l_br * tanxr*cosxr + l_bt * cosxr;  // ez a régi
	return     l_br * tanxr*cosxr + l_bt * cosxr;
}


// coordinate: R{CYL} -> R{TOR}
__device__ double cyl2tor_coord(double l_r, double l_t){
	double tanxr = l_t/l_r;

//	double cosxr = 1/sqrt(1+tanxr*tanxr);

	// radT = radC/cos(phi)
	return    l_r * sqrt(1+tanxr*tanxr);
//	return l_r;
}

/*
__device__ void cyl2tor(double l_br, double* l_bt, double l_r, double l_t, double* ans){

	double tanxr = l_t/l_r;
	double cosxr = 1/sqrt(1+tanxr*tanxr);


	// radT = radC * cos(phi) + torC * sin(phi)
	ans[0]  =  l_br * cosxr    + l_bt * tanxr*cosxr;

	// torT = -radC *   sin(phi)  + torC * cos(phi)
	ans[1]  = -l_br * tanxr*cosxr + l_br * cosxr;

	// radT = radC/cos(phi)
	ans[2]  = l_r/cosxr;

}*/


/*
            tan x
sin x  = --------------
	              1/2
          /       2   \
          \1 + tan x  /




               1
cos x  = --------------
	              1/2
          /       2   \
          \1 + tan x  /





sin x  = tan x * cos x

*/      
