__global__ void taiga(TaigaGlobals *g, TaigaCommons *c, double *service_var){
    // thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    service_var[6] = 3.1415926535897932456;
    
    //if (c->step_counter == 0)    g->detcellid[idx] = -1;
    
    if (g->detcellid[idx] == -1){
        //TaigaLocals *l;
        double generalised_coordinates[6];
        int detcellid = g->detcellid[idx];
        
        generalised_coordinates[0] = g->rad[idx];
        generalised_coordinates[1] = g->z[idx];
        generalised_coordinates[2] = g->tor[idx];
        generalised_coordinates[3] = g->vrad[idx];
        generalised_coordinates[4] = g->vz[idx];
        generalised_coordinates[5] = g->vtor[idx];
        
        g->detcellid[idx] = traj(c, generalised_coordinates, detcellid);
        
        g->rad[idx]  = generalised_coordinates[0];
        g->z[idx]    = generalised_coordinates[1];
        g->tor[idx]  = generalised_coordinates[2];
        g->vrad[idx] = generalised_coordinates[3];
        g->vz[idx]   = generalised_coordinates[4];
        g->vtor[idx] = generalised_coordinates[5];
    }
    service_var[0] = 42.24;
    service_var[3] =g->particle_number;
    service_var[4] =c->max_step_number;
    service_var[5] =c->step_counter;
    service_var[6]=g->rad[2];
}

__global__ void cuda_service_test(double *service_var){
    service_var[9] = 3.1415926535897932456;
}
