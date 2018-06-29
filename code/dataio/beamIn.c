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

	double *Rv, V_temp;

	int nrv = vectorReader0(&Rv, "input/manual_profile/rad.dat");
	printf("SIZEOF: %d\n",nrv);
	
	if(nrv<beam_number){
	    printf("ERROR: There is not enough beamlet in source\n");
    }

	double diam = (double)diameter / 1000;	
	double deflH = (double)deflH_degree/180*PI;
	double deflV = (double)deflV_degree/180*PI;	
    //printf("Angles:\t%lf;%lf\n",deflH,deflV);
    
	V_temp = sqrt(2*energy*1000*eperm);
	
	/* one-ion beam */ 
	if (beam_number == 1){	
		/* set position of particles */	
		XR[0]=Rv[0];
		XZ[0]=0.;
		XT[0]=0.;		

		/* set velocity of particles */
		
		VR[0] = -V_temp*cos(deflH)*cos(deflV);
		VZ[0] =  V_temp*sin(deflH);
		VT[0] =  V_temp*cos(deflH)*sin(deflV);
	/* more ion in beam */
	}else{
		/* initialize random generator */
		srand ( time(NULL) );
		for (i=0;i<beam_number;++i){
			/* set position of particles */
			XR[i]=Rv[i];
			do{
				XZ[i]=(double)(rand()-RAND_MAX/2)/RAND_MAX*diam;
				XT[i]=(double)(rand()-RAND_MAX/2)/RAND_MAX*diam;
			}while ((XZ[i]*XZ[i]+XT[i]*XT[i])>=(diam/2)*(diam/2));
			
			/* toroidal deflection */
			
			XT[i] += tan(deflV) * ($R_defl - XR[i]);
			
			/* set velocity of particles */
			VR[i] = -V_temp*cos(deflH)*cos(deflV);
			VZ[i] =  V_temp*sin(deflH);
			VT[i] =  V_temp*cos(deflH)*sin(deflV);		
			
			
		}
	}
	return beam_number;

}


 
