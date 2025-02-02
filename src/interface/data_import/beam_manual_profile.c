#include <stdlib.h>
#include <math.h>
#include "utils/physics.h"

void load_beam_manual(TaigaGlobals *g, BeamProp *beam, ShotProp *shot, RunProp *run){
    double *xr,*xz,*xt,*vr,*vz,*vt;
    long profile_length = -1;

    if (strcmp(run->ion_source_file, "")){
        if (!strcmp(run->io_coordinate_order, "rtz") || !strcmp(run->io_coordinate_order, "RTZ")){
            profile_length = read_matrix_column(&g->rad, run->ion_source_file, 1);
            read_matrix_column(&g->z, run->ion_source_file, 3);
            read_matrix_column(&g->tor, run->ion_source_file, 2);
            if(read_matrix_column(&g->vrad, run->ion_source_file, 4)){
                read_matrix_column(&g->vz, run->ion_source_file, 6);
                read_matrix_column(&g->vtor, run->ion_source_file, 5);
            }else{
                if (run->debug==1)  printf("Velocity is calculated from beam energy\n");
                double speed = calculate_speed(beam->energy, get_mass(beam->species, beam->charge));
                for (long i=0; i<profile_length; ++i){
                    g->vrad[i] = -speed*cos(beam->vertical_deflection)*cos(beam->toroidal_deflection);
                    g->vz[i]   =  speed*sin(beam->vertical_deflection);
                    g->vtor[i] =  speed*cos(beam->vertical_deflection)*sin(beam->toroidal_deflection);
                }
            }

        }else if (!strcmp(run->io_coordinate_order, "rzt") || !strcmp(run->io_coordinate_order, "RZT")){
            profile_length = read_matrix_column(&g->rad, run->ion_source_file, 1);
            read_matrix_column(&g->z, run->ion_source_file, 2);
            read_matrix_column(&g->tor, run->ion_source_file, 3);
            if (read_matrix_column(&g->vrad, run->ion_source_file, 4)){
                read_matrix_column(&g->vz, run->ion_source_file, 5);
                read_matrix_column(&g->vtor, run->ion_source_file, 6);
            }else{
                if (run->debug==1)  printf("Velocity is calculated from beam energy\n");
                double speed = calculate_speed(beam->energy, get_mass(beam->species, beam->charge));
                printf("energy: %lf keV,\tmass: %lf amu\tspeed: %lf m/s\n", beam->energy, get_mass(beam->species, beam->charge), speed);
                for (long i=0; i<profile_length; ++i){
                    g->vrad[i] = -speed*cos(beam->vertical_deflection)*cos(beam->toroidal_deflection);
                    g->vz[i]   =  speed*sin(beam->vertical_deflection);
                    g->vtor[i] =  speed*cos(beam->vertical_deflection)*sin(beam->toroidal_deflection);
                }
            }
            
        }else{
            printf("Invalid input format. Reading directly.");
            strcpy(run->ion_source_file, "");
        }
    }
    if (!strcmp(run->ion_source_file, "")){
        profile_length = read_vector(&g->rad, "input", "manual_profile", "rad.dat");
        read_vector(&g->z, "input", "manual_profile", "z.dat");
        read_vector(&g->tor, "input", "manual_profile", "tor.dat");
        read_vector(&g->vrad, "input", "manual_profile", "vrad.dat");
        read_vector(&g->vz, "input", "manual_profile", "vz.dat");
        read_vector(&g->vtor, "input", "manual_profile", "vtor.dat");
    }
    
    run->particle_number = profile_length;
    if (run->debug==1)  printf("Number of particles are modified to %d.\n", run->particle_number);
}
