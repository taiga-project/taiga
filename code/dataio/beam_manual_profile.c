#include <stdlib.h>
#include <math.h>

void load_beam(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, beam_prop beam, shot_prop shot, run_prop run){
    double *xr,*xz,*xt,*vr,*vz,*vt;
    int profile_length = -1;

    if (strcmp(run.ion_source_file, "")){
        if (!strcmp(run.io_coordinate_order, "rtz") || !strcmp(run.io_coordinate_order, "RTZ")){
            profile_length = read_matrix_column(&xr, run.ion_source_file, 1);
            read_matrix_column(&xz, run.ion_source_file, 3);
            read_matrix_column(&xt, run.ion_source_file, 2);
            read_matrix_column(&vr, run.ion_source_file, 4);
            read_matrix_column(&vz, run.ion_source_file, 6);
            read_matrix_column(&vt, run.ion_source_file, 5);

        }else if (!strcmp(run.io_coordinate_order, "rzt") || !strcmp(run.io_coordinate_order, "RZT")){
            profile_length = read_matrix_column(&xr, run.ion_source_file, 1);
            read_matrix_column(&xz, run.ion_source_file, 2);
            read_matrix_column(&xt, run.ion_source_file, 3);
            read_matrix_column(&vr, run.ion_source_file, 4);
            read_matrix_column(&vz, run.ion_source_file, 5);
            read_matrix_column(&vt, run.ion_source_file, 6);
            
        }else{
            printf("Invalid input format. Reading directly.");
            strcpy(run.ion_source_file, "");
        }
    }
    if (!strcmp(run.ion_source_file, "")){
        profile_length = read_vector(&xr, "input", "manual_profile", "rad.dat");
        read_vector(&xz, "input", "manual_profile", "z.dat");
        read_vector(&xt, "input", "manual_profile", "tor.dat");
        read_vector(&vr, "input", "manual_profile", "vrad.dat");
        read_vector(&vz, "input", "manual_profile", "vz.dat");
        read_vector(&vt, "input", "manual_profile", "vtor.dat");
    }
    
    run.particle_number = profile_length;
    if (run.debug==1)  printf("Number of particles are modified to %d.\n", run.particle_number);
    
    for (int i=0; i<run.particle_number; i++){
        XR[i]=xr[i];
        XZ[i]=xz[i];
        XT[i]=xt[i];
        VR[i]=vr[i];
        VZ[i]=vz[i];
        VT[i]=vt[i];
    }
}
