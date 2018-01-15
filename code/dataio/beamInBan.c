#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/*
diameter:	in mm
energy:		in keV
mass:		in AMU

*/

// set beam inline parameters
int beamIn(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, double energy, double eperm, int beam_number, char *shotname, double diameter, double deflH_degree, double deflV_degree){
	int i;

	double diam = (double)diameter / 1000;	
	double vmax = sqrt(2*energy*1000*eperm);

	/* one-ion beam */ 
	if (beam_number == 1){	
		/* set position of particles */	
		XR[0]=R_midions;
		XZ[0]=0.;
		XT[0]=0.;

		/* set velocity of particles */
		VR[0]=0.;
		VZ[0]=sqrt(2*energy*1000*eperm);
		VT[0]=0.;
	/* more ion in beam */
	}else{
		/* initialize random generator */
		srand ( time(NULL) );
		for (i=0;i<beam_number;++i){
			/* set position of particles */
			XR[i]=R_midions+(double)(rand()-RAND_MAX/2)/RAND_MAX*diam;
			XZ[i]=0;
			XT[i]=0;
			
			/* set velocity of particles */
			VR[i]=-0;
			VZ[i]=vmax*cos((double)(rand()-RAND_MAX/2)/RAND_MAX*PI);
			VT[i]=vmax*sin((double)(rand()-RAND_MAX/2)/RAND_MAX*PI);
		}
	}
	printf("Banana Joe\n");
	return beam_number;

}
