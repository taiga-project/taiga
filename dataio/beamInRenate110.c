#include <stdlib.h>
#include <math.h>



/*
diameter:	in mm
energy:		in keV
mass:		in AMU

*/
int get_array_size(double *array);
double linear_interpolate(double *x_vector, int x_length, double *y_vector, int y_length, double x_value);

// set beam inline parameters
void beamIn(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, double energy, double eperm, int beam_number, char *shotname, double diameter, double deflH_degree, double deflV_degree){
	int i;

	double *prof_r, *prof_d, Vabs, ionisation_yeald;

	int prof_r_length = vectorReader(&prof_r, "input/ionProf", shotname, "rad.dat");
	int prof_d_length = vectorReader(&prof_d, "input/ionProf", shotname, "ionyeald.dat");
	
	if (prof_r_length != prof_d_length){
		printf("ERROR: Length of PROF_R and RPOF_D are different!\n");
	}

	double diam = diameter / 1000.0;	
	double deflH = deflH_degree/180*PI;
	double deflV = deflV_degree/180*PI;	
    //printf("Angles:\t%lf;%lf\n",deflH,deflV);
    
	Vabs = sqrt(2*energy*1000*eperm);
	
	/* one-ion beam */ 

	/* initialize random generator */
	srand ( time(NULL) );
	for (i=0;i<beam_number;++i){
		/* set position of particles */
		do{
		    ionisation_yeald = (double)rand()/RAND_MAX;
		    XR[i] = linear_interpolate(prof_d, prof_d_length, prof_r, prof_r_length, ionisation_yeald);
	    }while (isnan(XR[i]));
		do{
			XZ[i]=(double)(rand()-RAND_MAX/2)*diam;
			XT[i]=(double)(rand()-RAND_MAX/2)*diam;
		    printf("(%d.) %lf %lf %lf\n",i, XR[i],XZ[i], XT[i]);
		}while ((XZ[i]*XZ[i]+XT[i]*XT[i])>=(diam/2)*(diam/2));
		
		/* toroidal deflection */
		
		XT[i] += tan(deflV) * ($R_defl - XR[i]);
		
		/* set velocity of particles */
		VR[i] = -Vabs*cos(deflH)*cos(deflV);
		VZ[i] =  Vabs*sin(deflH);
		VT[i] =  Vabs*cos(deflH)*sin(deflV);		
			
			
		
	}

}

double linear_interpolate(double *x_vector, int x_length, double *y_vector, int y_length, double x_value){

	int i;

/*	int x_length = get_array_size(x_vector);
	int y_length = get_array_size(y_vector);*/

	if (x_length != y_length)		printf("ERROR: in interpolation. Two input vectors have different length.");
	
	for (i=1; (i<x_length) && (x_vector[i-1]>x_value); i++);//{printf("%lf %lf |",x_vector[i-1],x_value);}
	
	if(i>1)--i;
//	printf("ii%ld",i);
	return y_vector[i] - (y_vector[i]-y_vector[i-1])*(x_value-x_vector[i-1])/(x_vector[i]-x_vector[i-1]);
	
}

int get_array_size(double *array){
	return (int)(sizeof(array)/sizeof(array[0]));
} 
