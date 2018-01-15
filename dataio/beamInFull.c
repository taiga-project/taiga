#include <stdlib.h>
#include <math.h>



/*
diameter:	in mm
energy:		in keV
mass:		in AMU

*/

// set beam inline parameters
int beamIn(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, double energy, double eperm, int beam_number, char *shotname, double diameter, double deflH_degree, double deflV_degree){

	int nrv = vectorReader0(&XR, "dataio/data/rad.dat");
	vectorReader0(&XZ, "dataio/data/z.dat");
	vectorReader0(&XT, "dataio/data/tor.dat");

	vectorReader0(&VR, "dataio/data/vrad.dat");
	vectorReader0(&VZ, "dataio/data/vz.dat");
	vectorReader0(&VT, "dataio/data/vtor.dat");

	printf("SIZEOF: %d\n",nrv);
	
	return nrv;

}
