__global__ void taiga(double timestep, int NR, int NZ, double eperm, double **spline_brad, double **spline_bz, double **spline_btor, double **spline_grid, double **position_all, double **speed_all, double *detector_geometry, int *detcellid, int N_step, double *service_var, int step_i,
                      int *prof_size, double *prof_r, double *prof_d, double *profx_r, double *profx_d){
    // thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (step_i == 0){
        detcellid[idx] = -1;
        if (prof_size[0] > 0){
            // GENERATE PROFILE
        }
    }
    
    if (detcellid[idx] == -1){
        double *spline_grid_rad, *spline_grid_z;
        spline_grid_rad = spline_grid[0];
        spline_grid_z   = spline_grid[1];

        double position[3], speed[3];

        position[0] = position_all[0][idx];
        position[1] = position_all[1][idx];
        position[2] = position_all[2][idx];

        speed[0] = speed_all[0][idx];
        speed[1] = speed_all[1][idx];
        speed[2] = speed_all[2][idx];
        
        detcellid[idx] = traj(timestep, spline_grid_rad, NR, spline_grid_z, NZ, position, speed, spline_brad, spline_bz, spline_btor, eperm, detector_geometry, N_step, detcellid[idx]);

        position_all[0][idx] = position[0];
        position_all[1][idx] = position[1];
        position_all[2][idx] = position[2];

        speed_all[0][idx] = speed[0];
        speed_all[1][idx] = speed[1];
        speed_all[2][idx] = speed[2];
    }

    service_var[0] = 42.24;
}

__global__ void taiga(double timestep, int NR, int NZ, double eperm, double **spline_brad, double **spline_bz, double **spline_btor, double **spline_erad, double **spline_ez, double **spline_etor, 
                      double **spline_grid, double **position_all, double **speed_all, double *detector_geometry, int *detcellid, int N_step, double *service_var, int step_i,
                      int *prof_size, double *prof_r, double *prof_d, double *profx_r, double *profx_d){
    // thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (step_i == 0){
        detcellid[idx] = -1;
        if (prof_size[0] > 0){
            // GENERATE PROFILE
        }
    }
    
    if (detcellid[idx] == -1){
        double *spline_grid_rad, *spline_grid_z;
        spline_grid_rad = spline_grid[0];
        spline_grid_z   = spline_grid[1];
        
        double position[3], speed[3];

        position[0] = position_all[0][idx];
        position[1] = position_all[1][idx];
        position[2] = position_all[2][idx];

        speed[0] = speed_all[0][idx];
        speed[1] = speed_all[1][idx];
        speed[2] = speed_all[2][idx];

        detcellid[idx] = traj(timestep, spline_grid_rad, NR, spline_grid_z, NZ, position, speed, spline_brad, spline_bz, spline_btor, spline_erad, spline_ez, spline_etor, eperm, detector_geometry, N_step, detcellid[idx]);

        position_all[0][idx] = position[0];
        position_all[1][idx] = position[1];
        position_all[2][idx] = position[2];

        speed_all[0][idx] = speed[0];
        speed_all[1][idx] = speed[1];
        speed_all[2][idx] = speed[2];
    }

    service_var[0] = 42.24;
}
