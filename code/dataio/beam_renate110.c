#include <stdlib.h>
#include <math.h>

// set beam inline parameters
void load_beam(double **X, beam_prop beam, shot_prop shot, run_prop run){
    int i, prof_size[2];
    double *prof_r, *prof_d, *profx_r, *profx_d, Vabs, ionisation_yeald, xsec_rad, xsec_ang;
    
    char* shotname = concat(shot.shotnumber, "_", shot.time);
    
    beam_profile prof;
    init_ion_profile(shotname, &prof);
    load_ion_profile(shotname, &prof);
    
    Vabs = sqrt(2 * beam.energy*1000*ELEMENTARY_CHARGE/ beam.mass/ AMU);
    
    /* cross section normalisation */
    /*if (prof_size[1] > 0){
        for (i=0; i<prof_size[1]; i++){
            profx_d[i] /= profx_r[i];
        }
    }*/
    
    /* initialize random generator */
    srand ( time(NULL) );
    
    for (i=0; i<run.particle_number; i++){
        /* set position of particles */
        do{
            ionisation_yeald = (double)rand()/RAND_MAX;
            X[0][i] = linear_interpolate(prof.radial.profile, prof.radial.N, prof.radial.grid, prof.radial.N, ionisation_yeald);
        }while (isnan(X[0][i])||X[0][i]<0);
        do{
            if (prof_size[1] <= 0){
                X[1][i]=(double)(rand()-RAND_MAX/2)/RAND_MAX*beam.diameter;
                X[2][i]=(double)(rand()-RAND_MAX/2)/RAND_MAX*beam.diameter;
            }else{
                ionisation_yeald = (double)rand()/RAND_MAX;
                xsec_ang = (double)rand()/RAND_MAX*2*PI;
                xsec_rad = linear_interpolate(prof.cross_section.profile, prof.cross_section.N, prof.cross_section.grid, prof.cross_section.N, ionisation_yeald)*(beam.diameter/2);
                X[1][i]= sin(xsec_ang) * xsec_rad;
                X[2][i]= cos(xsec_ang) * xsec_rad;
            }
        }while ((X[1][i]*X[1][i]+X[2][i]*X[2][i])>=(beam.diameter/2)*(beam.diameter/2));
        
        /* deflection */
        X[1][i] += tan(beam.vertical_deflection) * ($R_defl - X[0][i]);
        X[2][i] += tan(beam.toroidal_deflection) * ($R_defl - X[0][i]);
        
        /* set velocity of particles */
        X[3][i] = -Vabs*cos(beam.vertical_deflection)*cos(beam.toroidal_deflection);
        X[4][i] =  Vabs*sin(beam.vertical_deflection);
        X[5][i] =  Vabs*cos(beam.vertical_deflection)*sin(beam.toroidal_deflection);
    }
}

void init_ion_profile(char* shotname, beam_profile* prof){
    int prof_r_length = read_vector(&(prof->radial.grid),    "input/ionProf", shotname, "rad.dat");
    int prof_d_length = read_vector(&(prof->radial.profile), "input/ionProf", shotname, "ionyeald.dat");    
    
    int profx_r_length = read_vector(&(prof->cross_section.grid),    "input/ionProf", shotname, "xrad.dat", false);
    int profx_d_length = read_vector(&(prof->cross_section.profile), "input/ionProf", shotname, "xionyeald.dat", false);
    
    if (prof_r_length <= 1){
        printf("ERROR: Invalid length of PROF_R!\n");
        exit(0);
    }
    
    if (prof_r_length == prof_d_length){
        prof->radial.N = prof_r_length;
    }else{
        printf("ERROR: Length of PROF_R and PROF_D are different!\n");
        exit(0);
    }

    if (profx_r_length == profx_d_length){
        if (profx_r_length <= 1){
            printf("Cross section beam profile: OFF\n");            
            prof->cross_section.N = 0;
        }else{
            printf("Cross section beam profile: ON\n");
            prof->cross_section.N = profx_r_length;
        }
    }else{
        printf("WARNING: Length of PROFX_R and PROFX_D are different!\nCross section beam profile: OFF\n");
        prof->cross_section.N = 0;
    }    
}

void load_ion_profile(char* shotname, beam_profile* prof){

    double *local_prof_r, *local_prof_d, *local_profx_r, *local_profx_d;
    int i;

    read_vector(&local_prof_r, "input/ionProf", shotname, "rad.dat");
    read_vector(&local_prof_d, "input/ionProf", shotname, "ionyeald.dat");
    for (i=0; i<prof->radial.N; ++i){
        prof->radial.grid[i]    = local_prof_r[i];
        prof->radial.profile[i] = local_prof_d[i];
    }
    
    if (prof->cross_section.N > 1){
        read_vector(&local_profx_r, "input/ionProf", shotname, "xrad.dat", false);
        read_vector(&local_profx_d, "input/ionProf", shotname, "xionyeald.dat", false);
        
        for (i=0; i<prof->cross_section.N; i++){
            prof->cross_section.grid[i]    = local_profx_r[i];
            prof->cross_section.profile[i] = local_profx_d[i];
        }
    }
}
