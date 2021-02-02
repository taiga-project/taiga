__device__ double device_linear_interpolate(double *x_vector, int x_length, double *y_vector, int y_length, double x_value){
    int i;
    for (i=1; (i<x_length) && (x_vector[i-1]>x_value); ++i);
    if(i>1){--i;}else{i=1;}
    return y_vector[i] - (y_vector[i]-y_vector[i-1])*(x_value-x_vector[i-1])/(x_vector[i]-x_vector[i-1]);
}


/*(double beam.diameter, double beam.energy, double beam.vertical_deflection, double beam.toroidal_deflection,
                                double **position_all, double **speed_all, double eperm, int *prof_size, double *prof.radial.grid, double *prof.radial.profile, double *prof.cross_section.grid, double *prof.cross_section.profile)*/
__global__ void generate_coords(TaigaGlobals *g, TaigaCommons *s, BeamProp beam, BeamProfile prof){

    // thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    int i;
    double Vabs, ionisation_yeald, xsec_rad, xsec_ang;
    double XR, XZ, XT;
    
    Vabs = sqrt(2*beam.energy*1000*s->eperm);
    
    // cross section normalisation 
    /*if (prof.cross_section.N > 0){
        for (i=0; i<prof.cross_section.N; ++i){
            prof.cross_section.profile[i] /= prof.cross_section.grid[i];
        }
    }*/
    
    // initialize random generator 
    curandState state;
    curand_init((unsigned long long)clock() + idx, 0, 0, &state);
    
    // set position of particles 
    do{
        ionisation_yeald = curand_uniform_double(&state);
        XR = device_linear_interpolate(prof.radial.profile, prof.radial.N, prof.radial.grid, prof.radial.N, ionisation_yeald);
        g->rad[idx] = XR;
    }while (isnan(XR)||XR<0);
    do{
        //if (prof.cross_section.N <= 0){
            XZ = (curand_uniform_double(&state)-0.5)*beam.diameter;
            XT = (curand_uniform_double(&state)-0.5)*beam.diameter;
            g->z[idx] = XZ;
            g->tor[idx] = XT;
        //}else{
         //#   ionisation_yeald = curand_uniform_double(&state);
         //   xsec_ang = curand_uniform_double(&state)*2*PI;
         //   xsec_rad = linear_interpolate(prof.cross_section.profile, prof.cross_section.N, prof.cross_section.grid, prof.cross_section.N, ionisation_yeald)*(beam.diameter/2);
         //   XZ[i]= sin(xsec_ang) * xsec_rad;
         //   XT[i]= cos(xsec_ang) * xsec_rad;
        //}
    }while ((XZ*XZ+XT*XT)>=(beam.diameter/2)*(beam.diameter/2));
    
    
    // deflection 
    g->z[idx] += tan(beam.vertical_deflection) * (beam.deflection_coordinate - XR);
    g->tor[idx] += tan(beam.toroidal_deflection) * (beam.deflection_coordinate - XR);
    
    // set velocity of particles
    g->vrad[idx] = -Vabs*cos(beam.vertical_deflection)*cos(beam.toroidal_deflection);
    g->vz[idx] =  Vabs*sin(beam.vertical_deflection);
    g->vtor[idx] =  Vabs*cos(beam.vertical_deflection)*sin(beam.toroidal_deflection);
}
