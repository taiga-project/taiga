__global__ void taiga(TaigaGlobals *g, TaigaCommons *c, double *service_var){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;    // thread index
    
    if (g->detcellid[idx] == -1){
        double generalised_coordinates[6];
        int detcellid = g->detcellid[idx];
        
        generalised_coordinates[0] = g->rad[idx];
        generalised_coordinates[1] = g->z[idx];
        generalised_coordinates[2] = g->tor[idx];
        generalised_coordinates[3] = g->vrad[idx];
        generalised_coordinates[4] = g->vz[idx];
        generalised_coordinates[5] = g->vtor[idx];
        
        g->detcellid[idx] = calculate_trajectory(c, generalised_coordinates, detcellid);
        
        g->rad[idx]  = generalised_coordinates[0];
        g->z[idx]    = generalised_coordinates[1];
        g->tor[idx]  = generalised_coordinates[2];
        g->vrad[idx] = generalised_coordinates[3];
        g->vz[idx]   = generalised_coordinates[4];
        g->vtor[idx] = generalised_coordinates[5];
    }
    service_var[0] = 42.24;
}

__global__ void cuda_service_test(double *service_var){
    service_var[9] = 3.1415926535897932456;
}
