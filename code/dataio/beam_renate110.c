#include <stdlib.h>
#include <math.h>
#include beam.h

// set beam inline parameters
void load_beam(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, beam_prop beam, shot_prop shot, run_prop run){
    int i;
    double prof_size[2], *prof_r, *prof_d, *profx_r, *profx_d, Vabs, ionisation_yeald, xsec_rad, xsec_ang;
    
    char* shotname = concat(shot.shotnumber, "_", shot.time);

    load_ion_profile(shotname, prof_size, prof_r, prof_d, profx_r, profx_d);
    
    Vabs = sqrt(2 * beam.energy*1000*ELEMENTARY_CHARGE/ beam.mass/ AMU);
    
    /* cross section normalisation */
    if (prof_size[1] > 0){
        for (i=0; i<prof_size[1]; i++){
            profx_d[i] /= profx_r[i];
        }
    }
    
    /* initialize random generator */
    srand ( time(NULL) );
    
    for (i=0; i<run.particle_number; i++){
        /* set position of particles */
        do{
            ionisation_yeald = (double)rand()/RAND_MAX;
            XR[i] = linear_interpolate(prof_d, prof_size[0], prof_r, prof_size[0], ionisation_yeald);
        }while (isnan(XR[i])||XR[i]<0);
        do{
            if (prof_size[1] <= 0){
                XZ[i]=(double)(rand()-RAND_MAX/2)/RAND_MAX*beam.diameter;
                XT[i]=(double)(rand()-RAND_MAX/2)/RAND_MAX*beam.diameter;
            }else{
                ionisation_yeald = (double)rand()/RAND_MAX;
                xsec_ang = (double)rand()/RAND_MAX*2*PI;
                xsec_rad = linear_interpolate(profx_d, prof_size[1], profx_r, prof_size[1], ionisation_yeald)*(beam.diameter/2);
                XZ[i]= sin(xsec_ang) * xsec_rad;
                XT[i]= cos(xsec_ang) * xsec_rad;
            }
        }while ((XZ[i]*XZ[i]+XT[i]*XT[i])>=(beam.diameter/2)*(beam.diameter/2));
        
        /* toroidal deflection */
        XT[i] += tan(beam.toroidal_deflection) * ($R_defl - XR[i]);
        
        /* set velocity of particles */
        VR[i] = -Vabs*cos(beam.vertical_deflection)*cos(beam.toroidal_deflection);
        VZ[i] =  Vabs*sin(beam.vertical_deflection);
        VT[i] =  Vabs*cos(beam.vertical_deflection)*sin(beam.toroidal_deflection);
    }
}

void load_ion_profile(char* shotname, double *prof_size, double *prof_r, double *prof_d, double *profx_r, double *profx_d){
    int prof_r_length = read_vector(&prof_r, "input/ionProf", shotname, "rad.dat");
    int prof_d_length = read_vector(&prof_d, "input/ionProf", shotname, "ionyeald.dat");    
    
    int profx_r_length = read_vector(&profx_r, "input/ionProf", shotname, "xrad.dat", false);
    int profx_d_length = read_vector(&profx_d, "input/ionProf", shotname, "xionyeald.dat", false);
    
    if (prof_r_length == prof_d_length){
        prof_size[0] = prof_r_length;
    }else{
        printf("ERROR: Length of PROF_R and PROF_D are different!\n");
        prof_size[0] = -1;
    }

    if (profx_r_length == profx_d_length){
        prof_size[1] = profx_r_length;
    }else{
        printf("ERROR: Length of PROFX_R and PROFX_D are different!\n");
        prof_size[1] = -1;
    }

    if (prof_size[1] <= 0){
        printf("Cross section beam profile: OFF\n");
    }else{
        printf("Cross section beam profile: ON\n");
    }
    
}
