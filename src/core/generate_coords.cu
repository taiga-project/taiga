#include "maths.cuh"

__global__ void generate_coords(TaigaGlobals *globals, BeamProp beam, BeamProfile *prof){

    // thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    double speed, ionisation_yeald;
    double XR, XZ, XT;
    
    speed = sqrt(2 * beam.energy*1000*ELEMENTARY_CHARGE/ AMU/ beam.mass);
    
    // cross section normalisation 
//    if (prof->cross_length > 0){
//        for (i=0; i<prof->cross_length; ++i){
//            prof->cross_profile[i] /= prof->cross_grid[i];
//        }
//    }
    
    // initialize random generator 
    curandState state;
    curand_init((unsigned long long)clock() + idx, 0, 0, &state);
    
    // set position of particles 
    do{
        ionisation_yeald = curand_uniform_double(&state);
        XR = interpolate_from_vector(prof->radial_profile, prof->radial_grid, prof->radial_length, ionisation_yeald);
        globals->rad[idx] = XR;
    }while (isnan(XR)||XR<0);
    do{
        //if (prof->cross_length <= 0){
            XZ = (curand_uniform_double(&state)-0.5)*beam.diameter;
            XT = (curand_uniform_double(&state)-0.5)*beam.diameter;
            globals->z[idx] = XZ;
            globals->tor[idx] = XT;
        //}else{
        // double xsec_rad, xsec_ang;
         //#   ionisation_yeald = curand_uniform_double(&state);
         //   xsec_ang = curand_uniform_double(&state)*2*PI;
         //   xsec_rad = linear_interpolate(prof->cross_profile, prof->cross_length, prof->cross_grid, prof->cross_length, ionisation_yeald)*(beam.diameter/2);
         //   XZ[i]= sin(xsec_ang) * xsec_rad;
         //   XT[i]= cos(xsec_ang) * xsec_rad;
        //}
    }while ((XZ*XZ+XT*XT)>=(beam.diameter/2)*(beam.diameter/2));
    
    // deflection 
    globals->z[idx] += tan(beam.vertical_deflection) * (beam.deflection_coordinate - XR);
    globals->tor[idx] += tan(beam.toroidal_deflection) * (beam.deflection_coordinate - XR);
    
    // set velocity of particles*/
    globals->vrad[idx] = -speed*cos(beam.vertical_deflection)*cos(beam.toroidal_deflection);
    globals->vz[idx] =  speed*sin(beam.vertical_deflection);
    globals->vtor[idx] =  speed*cos(beam.vertical_deflection)*sin(beam.toroidal_deflection);
    globals->detcellid[idx] = CALCULATION_NOT_FINISHED;
}
