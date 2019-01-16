#include <stdlib.h>
#include <math.h>

// set beam inline parameters
int beamIn(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, double energy, double eperm, int beam_number, char *shotname, double diameter, double deflH_degree, double deflV_degree){
	int i;
	double *xr,*xz,*xt,*vr,*vz,*vt;
    vectorReader0(&xr, "input/manual_profile/rad.dat");
	vectorReader0(&xz, "input/manual_profile/z.dat");
	vectorReader0(&xt, "input/manual_profile/tor.dat");
	vectorReader0(&vr, "input/manual_profile/vrad.dat");
	vectorReader0(&vz, "input/manual_profile/vz.dat");
	vectorReader0(&vt, "input/manual_profile/vtor.dat");

	for (i=0;i<beam_number;++i){
		/* set position of particles */
		XR[i]=xr[i];
		XZ[i]=xz[i];
		XT[i]=xt[i];
		
		/* set velocity of particles */
		VR[i]=vr[i];
		VZ[i]=vz[i];
		VT[i]=vt[i];
	}
	return beam_number;
}
