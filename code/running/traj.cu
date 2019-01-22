__device__ /*double int*/ void copy_local_field(double *r_grid, int NR, double *z_grid, int NZ, double position_rad, double position_z, int *local_spline_indices, double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,  double **br_ptr, double **bz_ptr, double **bt_ptr){
	int rci, zci;
	int i, i2;
	
	for(rci=0;(r_grid[rci+1]<position_rad)&&(rci<NR-1);rci++){;}
	
	for(zci=0;(z_grid[zci+1]<position_z)&&(zci<NR-1);zci++){;}
	
	// Particle leave out the cell
	if ((local_spline_indices[0]!=rci)||(local_spline_indices[1]!=zci)){
		local_spline_indices[0]=rci;
		local_spline_indices[1]=zci;
	
		for(i=0;i<16;i++){
	
			i2 = (local_spline_indices[0])*(NZ-1)+local_spline_indices[1];
		
			local_spline_brad[i]=br_ptr[i][i2];
			local_spline_bz[i]=bz_ptr[i][i2];
			local_spline_btor[i]=bt_ptr[i][i2];
		
		}
	}
	
	//return i2;	
}

__device__ /*double int*/ void copy_local_field(double *r_grid, int NR, double *z_grid, int NZ, double position_rad, double position_z, int *local_spline_indices, double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,  double **br_ptr, double **bz_ptr, double **bt_ptr,
										double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,  double **er_ptr, double **ez_ptr, double **et_ptr){
	int rci, zci;
	int i, i2;
	
	for(rci=0;(r_grid[rci+1]<position_rad)&&(rci<NR-1);rci++){;}
	
	for(zci=0;(z_grid[zci+1]<position_z)&&(zci<NR-1);zci++){;}
	
	// Particle leave out the cell
	if ((local_spline_indices[0]!=rci)||(local_spline_indices[1]!=zci)){
		local_spline_indices[0]=rci;
		local_spline_indices[1]=zci;
	
		for(i=0;i<16;i++){
	
			i2 = (local_spline_indices[0])*(NZ-1)+local_spline_indices[1];
		
			local_spline_brad[i]=br_ptr[i][i2];
			local_spline_bz[i]=bz_ptr[i][i2];
			local_spline_btor[i]=bt_ptr[i][i2];
			local_spline_erad[i]=er_ptr[i][i2];
			local_spline_ez[i]=ez_ptr[i][i2];
			local_spline_etor[i]=et_ptr[i][i2];
		
		}
	}
	
	//return i2;	
}

__device__ double calculate_local_field(double *local_spline, double dr, double dz){

	/* MATLAB CODE:
	sample2(3) =c11(bs1,bs2)*dsx^3*dsy^3 + c12(bs1,bs2)*dsx^3*dsy^2 + c13(bs1,bs2)*dsx^3*dsy + c14(bs1,bs2)*dsx^3 + ...
				c21(bs1,bs2)*dsx^2*dsy^3 + c22(bs1,bs2)*dsx^3*dsy^2 + c23(bs1,bs2)*dsx^2*dsy + c24(bs1,bs2)*dsx^2 + ...
				c31(bs1,bs2)*dsx  *dsy^3 + c32(bs1,bs2)*dsx  *dsy^2 + c33(bs1,bs2)*dsx  *dsy + c34(bs1,bs2)*dsx    + ...
				c41(bs1,bs2)      *dsy^3 + c42(bs1,bs2)      *dsy^2 + c43(bs1,bs2)      *dsy + c44(bs1,bs2);*/

	double local_field = 0.0, tmp[16] ;
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			tmp[i*4+j] = local_spline[i*4+j]*pow(dr,3-i)*pow(dz,3-j);
		}
	}   
    
	for(int i=0;i<4;i++){
    	for(int j=0;j<4;j++){
			local_field+=tmp[i*4+j];
		}
	}   

	return local_field;
}


__device__ int traj(double *r_grid, int NR, double *z_grid, int NZ, double *position, double *speed, double **br_ptr, double **bz_ptr, double **bt_ptr, double eperm, double *detector_geometry, int N_step, int local_detcellid){

	// next grid
	int local_spline_indices[2];
	local_spline_indices[0]=-1;
	local_spline_indices[1]=-1;
		
	double local_spline_brad[16];
	double local_spline_bz[16];
	double local_spline_btor[16];
	
	double local_brad=0,local_bz,local_btor;
	double dr,dz;
	double position_rad_torus;
	
	double speed_rad, speed_z, speed_tor, speed_rad_prev, speed_z_prev, speed_tor_prev;
	double position_rad,  position_z,   position_tor, position_rad_prev,  position_z_prev,  position_tor_prev;

	double X[6];
	
	int finished = local_detcellid + 1;
	
	position_rad  = position[0];
	position_z  = position[1];
	position_tor  = position[2];
	
	speed_rad = speed[0];
	speed_z = speed[1];
	speed_tor = speed[2];
	
	int loopi;
	for (loopi=0;(loopi<N_step && (!finished));loopi++){
		// Get local magnetic field

		position_rad_torus = cyl2tor_coord(position_rad, position_tor);
		copy_local_field(r_grid,NR,z_grid,NZ,position_rad_torus,position_z,local_spline_indices,local_spline_brad,local_spline_bz,local_spline_btor,br_ptr,bz_ptr,bt_ptr);
		
		dr = position_rad_torus-r_grid[local_spline_indices[0]];
		dz = position_z-z_grid[local_spline_indices[1]];
	
		local_brad =  calculate_local_field(local_spline_brad,dr,dz);
		local_bz =  calculate_local_field(local_spline_bz,dr,dz);
		local_btor =  calculate_local_field(local_spline_btor,dr,dz);

		local_brad = cyl2tor_rad(local_brad, local_btor, position_rad, position_tor);
		local_btor = cyl2tor_field(local_brad, local_btor, position_rad, position_tor);
	
		// archivate coordinates
		position_rad_prev  = position_rad;	position_z_prev  = position_z;	position_tor_prev  = position_tor;
		speed_rad_prev = speed_rad;	speed_z_prev = speed_z;	speed_tor_prev = speed_tor;

		// new coordinates
		X[0] = position_rad;
		X[1] = position_z;
		X[2] = position_tor;

		// new speed components
		X[3] = speed_rad;
		X[4] = speed_z;
		X[5] = speed_tor;
	
		solve_diffeq(X, local_brad, local_bz, local_btor, eperm);
		
		// new coordinates
		position_rad = X[0];
		position_z   = X[1];
		position_tor = X[2];

		// new speed components
		speed_rad = X[3];
		speed_z   = X[4];
		speed_tor = X[5];
	
		// finished? (interpolation)
		finished = calculate_detection_position(position_rad, position_z, position_tor, position_rad_prev, position_z_prev, position_tor_prev, detector_geometry, position, speed_rad);
	}
	
	speed[0] = speed_rad;
	speed[1] = speed_z;
	speed[2] = speed_tor;

	if (!finished){
		position[0] = position_rad;
		position[1] = position_z;
		position[2] = position_tor;
	}
	if (finished){
		local_detcellid = 0;
	}
	return local_detcellid;
}

__device__ int traj(double *r_grid, int NR, double *z_grid, int NZ, double *position, double *speed, double **br_ptr, double **bz_ptr, double **bt_ptr, double **er_ptr, double **ez_ptr, double **et_ptr, double eperm, double *detector_geometry, int N_step, int local_detcellid){

	// next grid
	int local_spline_indices[2];
	local_spline_indices[0]=-1;
	local_spline_indices[1]=-1;
    
	double local_spline_brad[16];
	double local_spline_bz[16];
	double local_spline_btor[16];

	double local_spline_erad[16];
	double local_spline_ez[16];
	double local_spline_etor[16];
	
	double local_brad=0,local_bz=0,local_btor=0;
	double local_erad=0,local_ez=0,local_etor=0;
	double dr, dz;
	double position_rad_torus;
	
	double speed_rad, speed_z, speed_tor, speed_rad_prev, speed_z_prev, speed_tor_prev;
	double position_rad,  position_z,   position_tor, position_rad_prev,  position_z_prev,  position_tor_prev;

	double X[6];
	
	int finished = local_detcellid + 1;

	position_rad  = position[0];
	position_z    = position[1];
	position_tor  = position[2];
	
	speed_rad = speed[0];
	speed_z   = speed[1];
	speed_tor = speed[2];
	
	int loopi;
	for (loopi=0;(loopi<N_step && (!finished));loopi++){
		// Get local magnetic field

		position_rad_torus = cyl2tor_coord(position_rad, position_tor);
		copy_local_field(r_grid,NR,z_grid,NZ,position_rad_torus,position_z,local_spline_indices,local_spline_brad,local_spline_bz,local_spline_btor,br_ptr,bz_ptr,bt_ptr,local_spline_erad,local_spline_ez,local_spline_etor,er_ptr,ez_ptr,et_ptr);
		
		dr = position_rad_torus-r_grid[local_spline_indices[0]];
		dz = position_z-z_grid[local_spline_indices[1]];
	
		local_brad =  calculate_local_field(local_spline_brad,dr,dz);
		local_bz =  calculate_local_field(local_spline_bz,dr,dz);
		local_btor =  calculate_local_field(local_spline_btor,dr,dz);
		local_brad = cyl2tor_rad(local_brad, local_btor, position_rad, position_tor);
		local_btor = cyl2tor_field(local_brad, local_btor, position_rad, position_tor);
		
		local_erad =  calculate_local_field(local_spline_erad,dr,dz);
		local_ez =  calculate_local_field(local_spline_ez,dr,dz);
		local_etor =  calculate_local_field(local_spline_etor,dr,dz);
		local_erad = cyl2tor_rad(local_erad, local_etor, position_rad, position_tor);
		local_etor = cyl2tor_field(local_erad, local_etor, position_rad, position_tor);

		// archivate coordinates
		position_rad_prev  = position_rad;	position_z_prev  = position_z;	position_tor_prev  = position_tor;
		speed_rad_prev = speed_rad;	speed_z_prev = speed_z;	speed_tor_prev = speed_tor;

		// new coordinates
		X[0] = position_rad;
		X[1] = position_z;
		X[2] = position_tor;

		// new speed components
		X[3] = speed_rad;
		X[4] = speed_z;
		X[5] = speed_tor;
	
		solve_diffeq(X, local_brad, local_bz, local_btor, local_erad, local_ez, local_etor, eperm);
		
		// new coordinates
		position_rad = X[0];
		position_z   = X[1];
		position_tor = X[2];

		// new speed components
		speed_rad = X[3];
		speed_z   = X[4];
		speed_tor = X[5];
	
		// finished? (interpolation)
		finished = calculate_detection_position(position_rad, position_z, position_tor, position_rad_prev, position_z_prev, position_tor_prev, detector_geometry, position, speed_rad);
	}
	
	speed[0] = speed_rad;
	speed[1] = speed_z;
	speed[2] = speed_tor;
	
	if (!finished){
		position[0] = position_rad;
		position[1] = position_z;
		position[2] = position_tor;
	}

	if (finished){
		local_detcellid = 0;
	}

	return local_detcellid;
}
