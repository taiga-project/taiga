__global__ void taiga(device_global g, device_shared s, double *service_var){
//(double timestep, int NR, int NZ, double eperm, double **spline_brad, double **spline_bz, double **spline_btor, double **spline_grid, double **position_all, double **speed_all, double *detector_geometry, int *detcellid, int N_step, double *service_var, int step_i){
    // thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (s.step_counter == 0)    g.detcellid[idx] = -1;
    
    if (g.detcellid[idx] == -1){

        device_local l;

        for (int i=0; i<6; ++i){
            l.coords[i] = g.coords[idx][i];
        }

        l.detcellid = g.detcellid[idx];
        //g.detcellid[idx] = traj(l, s);
    }
    service_var[0] = 42.24;
}
