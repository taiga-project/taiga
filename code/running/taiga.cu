__global__ void taiga(/*device_global g, device_shared s,*/ double *service_var){
//(double timestep, int NR, int NZ, double eperm, double **spline_brad, double **spline_bz, double **spline_btor, double **spline_grid, double **position_all, double **speed_all, double *detector_geometry, int *detcellid, int N_step, double *service_var, int step_i){
    // thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    service_var[5] = 3.1415926535897932456;
    
    /*#if (s.step_counter == 0)    g.detcellid[idx] = -1;
    
    if (g.detcellid[idx] == -1){
        device_local l;
        //cudaMalloc((void **) &(l.coords),  6*sizeof(double));
        l.coords[0] = g.rad[idx];
        l.coords[1] = g.z[idx];
        l.coords[2] = g.tor[idx];
        l.coords[3] = g.vrad[idx];
        l.coords[4] = g.vz[idx];
        l.coords[5] = g.vtor[idx];
        
        l.detcellid = g.detcellid[idx];
        //#g.detcellid[idx] = traj(l, s);
        
        g.rad[idx]  = l.coords[0];
        g.z[idx]    = l.coords[1];
        g.tor[idx]  = l.coords[2];
        g.vrad[idx] = l.coords[3];
        g.vz[idx]   = l.coords[4];
        g.vtor[idx] = l.coords[5];
    }*/
    service_var[0] = 42.24;
}

__global__ void cuda_service_test(/*device_global g, device_shared s,*/ double *service_var){
    service_var[9] = 3.1415926535897932456;
}
