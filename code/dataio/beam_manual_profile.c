#include <stdlib.h>
#include <math.h>

void load_beam(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, beam_prop beam, shot_prop shot, run_prop run){
    double *xr,*xz,*xt,*vr,*vz,*vt;
    int profile_length = read_vector(&xr, "input", "manual_profile", "rad.dat");
    read_vector(&xz, "input", "manual_profile", "z.dat");
    read_vector(&xt, "input", "manual_profile", "tor.dat");
    read_vector(&vr, "input", "manual_profile", "vrad.dat");
    read_vector(&vz, "input", "manual_profile", "vz.dat");
    read_vector(&vt, "input", "manual_profile", "vtor.dat");
    
    run.particle_number = profile_length;
    if (shot.debug==1)  printf("Number of particles are modified to %d.\n", run.particle_number);
    
    for (int i=0; i<run.particle_number; i++){
        XR[i]=xr[i];
        XZ[i]=xz[i];
        XT[i]=xt[i];
        VR[i]=vr[i];
        VZ[i]=vz[i];
        VT[i]=vt[i];
    }
}
