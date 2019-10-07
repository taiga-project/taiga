__device__ double device_linear_interpolate(double *x_vector, int x_length, double *y_vector, int y_length, double x_value){
    int i;       
    for (i=1; (i<x_length) && (x_vector[i-1]>x_value); i++);    
    if(i>1){--i;}else{i=1;}    
    return y_vector[i] - (y_vector[i]-y_vector[i-1])*(x_value-x_vector[i-1])/(x_vector[i]-x_vector[i-1]);
}


/*(double beam_diameter, double beam_energy, double beam_vertical_deflection, double beam_toroidal_deflection,
                                double **position_all, double **speed_all, double eperm, int *prof_size, double *prof_r, double *prof_d, double *profx_r, double *profx_d)*/
__global__ void generate_coords(device_global &g, beam_prop b, beam_profile p){

    // thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    /*
    int i;
    double Vabs, ionisation_yeald, xsec_rad, xsec_ang;
    double XR, XZ, XT;
    
    Vabs = sqrt(2*beam_energy*1000*eperm);
    
    // cross section normalisation 
    if (prof_size[1] > 0){
        for (i=0; i<prof_size[1]; i++){
            profx_d[i] /= profx_r[i];
        }
    }
    
    // initialize random generator 
    curandState state;
	curand_init((unsigned long long)clock() + idx, 0, 0, &state);	
	    
    // set position of particles 
    do{
        ionisation_yeald = curand_uniform_double(&state);
        XR = device_linear_interpolate(prof_d, prof_size[0], prof_r, prof_size[0], ionisation_yeald);
        g->coords[0][idx] = XR;
    }while (isnan(XR)||XR<0);
    do{
        //if (prof_size[1] <= 0){
            XZ = (curand_uniform_double(&state)-0.5)*beam_diameter;
            XT = (curand_uniform_double(&state)-0.5)*beam_diameter;
            g->coords[1][idx] = XZ;
            g->coords[2][idx] = XT;
        }else{
         //#   ionisation_yeald = curand_uniform_double(&state);
         //   xsec_ang = curand_uniform_double(&state)*2*PI;
         //   xsec_rad = linear_interpolate(profx_d, prof_size[1], profx_r, prof_size[1], ionisation_yeald)*(beam_diameter/2);
         //   XZ[i]= sin(xsec_ang) * xsec_rad;
         //   XT[i]= cos(xsec_ang) * xsec_rad;
        }
    }while ((XZ*XZ+XT*XT)>=(beam_diameter/2)*(beam_diameter/2));
    
    
    // deflection 
    g->coords[1][idx] += tan(beam_vertical_deflection) * ($R_defl - XR);
    g->coords[2][idx] += tan(beam_toroidal_deflection) * ($R_defl - XR);
    
    // set velocity of particles
    g->coords[3][idx] = -Vabs*cos(beam_vertical_deflection)*cos(beam_toroidal_deflection);
    g->coords[4][idx] =  Vabs*sin(beam_vertical_deflection);
    g->coords[5][idx] =  Vabs*cos(beam_vertical_deflection)*sin(beam_toroidal_deflection);*/
}
