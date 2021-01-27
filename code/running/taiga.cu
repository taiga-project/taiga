__global__ void taiga(TaigaGlobals *g, TaigaCommons *s, double *service_var){
    // thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    service_var[6] = 3.1415926535897932456;
    
    //if (s->step_counter == 0)    g->detcellid[idx] = -1;
    
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

        //g->detcellid[idx] = traj(s, generalised_coordinates, detcellid);
        
        if (idx==10){
        /*service_var[6] = s->brad[0][0];
        service_var[7] = (double)s->grid_size[0]+(double)s->grid_size[1]/1000;
        service_var[3] = s->spline_rgrid[0];*//*
        service_var[4] = s->spline_zgrid[0];
        service_var[5] = s->spline_zgrid[1];*/
        }
        /*g->rad[idx]  = generalised_coordinates[0];
        g->z[idx]    = generalised_coordinates[1];
        g->tor[idx]  = generalised_coordinates[2];
        g->vrad[idx] = generalised_coordinates[3];
        g->vz[idx]   = generalised_coordinates[4];
        g->vtor[idx] = generalised_coordinates[5];*/
    }
    service_var[0] = 42.24;
    if (idx==10) service_var[2] =g->rad[2];
}

__global__ void cuda_service_test(double *service_var){
    service_var[9] = 3.1415926535897932456;
}
