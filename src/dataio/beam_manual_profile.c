#include <stdlib.h>
#include <math.h>

void load_beam(TaigaGlobals *g, BeamProp *beam, ShotProp *shot, RunProp *run){
    double *xr,*xz,*xt,*vr,*vz,*vt;
    int profile_length = -1;

    if (strcmp(run->ion_source_file, "")){
        if (!strcmp(run->io_coordinate_order, "rtz") || !strcmp(run->io_coordinate_order, "RTZ")){
            profile_length = read_matrix_column(&g->rad, run->ion_source_file, 1);
            read_matrix_column(&g->z, run->ion_source_file, 3);
            read_matrix_column(&g->tor, run->ion_source_file, 2);
            read_matrix_column(&g->vrad, run->ion_source_file, 4);
            read_matrix_column(&g->vz, run->ion_source_file, 6);
            read_matrix_column(&g->vtor, run->ion_source_file, 5);

        }else if (!strcmp(run->io_coordinate_order, "rzt") || !strcmp(run->io_coordinate_order, "RZT")){
            profile_length = read_matrix_column(&g->rad, run->ion_source_file, 1);
            read_matrix_column(&g->z, run->ion_source_file, 2);
            read_matrix_column(&g->tor, run->ion_source_file, 3);
            read_matrix_column(&g->vrad, run->ion_source_file, 4);
            read_matrix_column(&g->vz, run->ion_source_file, 5);
            read_matrix_column(&g->vtor, run->ion_source_file, 6);
            
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

void init_ion_profile(char* shotname, BeamProfile* prof){
    printf("Error: manual beam profile cannot call init_ion_profile");
    exit(1);
}
