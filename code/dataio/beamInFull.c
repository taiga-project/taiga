#include <stdlib.h>
#include <math.h>



/*
diameter:	in mm
energy:		in keV
mass:		in AMU

*/

// set beam inline parameters
int beamIn(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, double energy, double eperm, int beam_number, char *shotname, double diameter, double deflH_degree, double deflV_degree){

	int nrv = vectorReader0(&XR, "input/manual_profile/rad.dat");
	vectorReader0(&XZ, "input/manual_profile/z.dat");
	vectorReader0(&XT, "input/manual_profile/tor.dat");

	vectorReader0(&VR, "input/manual_profile/vrad.dat");
	vectorReader0(&VZ, "input/manual_profile/vz.dat");
	vectorReader0(&VT, "input/manual_profile/vtor.dat");

	printf("SIZEOF: %d\n",nrv);
	
	return nrv;

}
