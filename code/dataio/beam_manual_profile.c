#include <stdlib.h>
#include <math.h>

// set beam inline parameters
int beamIn(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, double energy, double eperm, int beam_number, char *shotname, double diameter, double deflH_degree, double deflV_degree){
	double *xr,*xz,*xt,*vr,*vz,*vt;
	int profile_length = read_vector(&xr, "input", "manual_profile", "rad.dat");
	read_vector(&xz, "input", "manual_profile", "z.dat");
	read_vector(&xt, "input", "manual_profile", "tor.dat");
	read_vector(&vr, "input", "manual_profile", "vrad.dat");
	read_vector(&vz, "input", "manual_profile", "vz.dat");
	read_vector(&vt, "input", "manual_profile", "vtor.dat");
	
	for (int i=0; i<beam_number; ++i){
		/* set position of particles */
		XR[i]=xr[i];
		XZ[i]=xz[i];
		XT[i]=xt[i];
		
		/* set velocity of particles */
		VR[i]=vr[i];
		VZ[i]=vz[i];
		VT[i]=vt[i];
	}
	if (profile_length > beam_number){
		beam_number = profile_length;
		printf("Number of particles are modified to %d.\n", beam_number);
	}
	return beam_number;
}
