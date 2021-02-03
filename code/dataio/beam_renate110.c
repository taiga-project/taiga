#include <stdlib.h>
#include <math.h>

// set beam inline parameters
void load_beam(TaigaGlobals *g, BeamProp beam, ShotProp shot, RunProp run){
    int i, prof_size[2];
    double *radial_grid, *radial_profile, *cross_section_grid, *cross_section_profile, speed, ionisation_yeald, xsec_rad, xsec_ang;
    
    char* shotname = concat(shot.shotnumber, "_", shot.time);
    
    BeamProfile prof;
    init_ion_profile(shotname, &prof);
    load_ion_profile(shotname, &prof);
    
    speed = sqrt(2 * beam.energy*1000*ELEMENTARY_CHARGE/ beam.mass/ AMU);
    
    /* cross section normalisation */
    /*if (prof_size[1] > 0){
        for (i=0; i<prof_size[1]; ++i){
            cross_section_profile[i] /= cross_section_grid[i];
        }
    }*/
    
    /* initialize random generator */
    srand ( time(NULL) );
    
    for (i=0; i<run.particle_number; ++i){
        /* set position of particles */
        do{
            ionisation_yeald = (double)rand()/RAND_MAX;
            g->rad[i] = linear_interpolate(prof.radial.profile, prof.radial.N, prof.radial.grid, prof.radial.N, ionisation_yeald);
        }while (isnan(g->rad[i])||g->rad[i]<0);
        do{
            if (prof_size[1] <= 0){
                g->z[i]   = (double)(rand()-RAND_MAX/2)/RAND_MAX*beam.diameter;
                g->tor[i] = (double)(rand()-RAND_MAX/2)/RAND_MAX*beam.diameter;
            }else{
                ionisation_yeald = (double)rand()/RAND_MAX;
                xsec_ang = (double)rand()/RAND_MAX*2*PI;
                xsec_rad = linear_interpolate(prof.cross_section.profile, prof.cross_section.N, prof.cross_section.grid, prof.cross_section.N, ionisation_yeald)*(beam.diameter/2);
                g->z[i]   = sin(xsec_ang) * xsec_rad;
                g->tor[i] = cos(xsec_ang) * xsec_rad;
            }
        }while ((g->z[i]*g->z[i]+g->tor[i]*g->tor[i])>=(beam.diameter/2)*(beam.diameter/2));
        
        /* deflection */
        g->z[i]   += tan(beam.vertical_deflection) * (beam.deflection_coordinate - g->rad[i]);
        g->tor[i] += tan(beam.toroidal_deflection) * (beam.deflection_coordinate - g->rad[i]);
        
        /* set velocity of particles */
        g->vrad[i] = -speed*cos(beam.vertical_deflection)*cos(beam.toroidal_deflection);
        g->vz[i]   =  speed*sin(beam.vertical_deflection);
        g->vtor[i] =  speed*cos(beam.vertical_deflection)*sin(beam.toroidal_deflection);
    }
    printf("Renate Particle number %d\n",run.particle_number);
}

void init_ion_profile(char* shotname, BeamProfile* prof){
    int radial_grid_length = read_vector(&(prof->radial.grid),    "input/ionProf", shotname, "rad.dat");
    int radial_profile_length = read_vector(&(prof->radial.profile), "input/ionProf", shotname, "ionyeald.dat");    
    
    int cross_section_grid_length = read_vector(&(prof->cross_section.grid),    "input/ionProf", shotname, "xrad.dat", false);
    int cross_section_profile_length = read_vector(&(prof->cross_section.profile), "input/ionProf", shotname, "xionyeald.dat", false);
    
    if (radial_grid_length <= 1){
        printf("ERROR: Invalid length of radial_grid!\n");
        exit(1);
    }
    
    if (radial_grid_length == radial_profile_length){
        prof->radial.N = radial_grid_length;
    }else{
        printf("ERROR: Length of radial_grid and radial_profile are different!\n");
        exit(1);
    }

    if (cross_section_grid_length == cross_section_profile_length){
        if (cross_section_grid_length <= 1){
            printf("Cross section beam profile: OFF\n");            
            prof->cross_section.N = 0;
        }else{
            printf("Cross section beam profile: ON\n");
            prof->cross_section.N = cross_section_grid_length;
        }
    }else{
        printf("WARNING: Length of cross_section_grid and cross_section_profile are different!\nCross section beam profile: OFF\n");
        prof->cross_section.N = 0;
    }    
}

void load_ion_profile(char* shotname, BeamProfile* prof){

    double *local_radial_grid, *local_radial_profile, *local_cross_section_grid, *local_cross_section_profile;
    int i;

    read_vector(&local_radial_grid, "input/ionProf", shotname, "rad.dat");
    read_vector(&local_radial_profile, "input/ionProf", shotname, "ionyeald.dat");
    for (i=0; i<prof->radial.N; ++i){
        prof->radial.grid[i]    = local_radial_grid[i];
        prof->radial.profile[i] = local_radial_profile[i];
    }
    
    if (prof->cross_section.N > 1){
        read_vector(&local_cross_section_grid, "input/ionProf", shotname, "xrad.dat", false);
        read_vector(&local_cross_section_profile, "input/ionProf", shotname, "xionyeald.dat", false);
        
        for (i=0; i<prof->cross_section.N; ++i){
            prof->cross_section.grid[i]    = local_cross_section_grid[i];
            prof->cross_section.profile[i] = local_cross_section_profile[i];
        }
    }
}
