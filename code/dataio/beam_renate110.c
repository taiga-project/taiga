#include <stdlib.h>
#include <math.h>

int get_array_size(double *array);
double linear_interpolate(double *x_vector, int x_length, double *y_vector, int y_length, double x_value);

// set beam inline parameters
void load_beam(double *XR, double *XZ, double *XT, double *VR, double *VZ, double *VT, beam_prop beam, shot_prop shot){
    int i;
    double *prof_r, *prof_d, *profx_r, *profx_d, Vabs, ionisation_yeald, xsec_rad, xsec_ang;
    
    char* shotname = concat(shot.shotnumber, "_", shot.time);

    int prof_r_length = read_vector(&prof_r, "input/ionProf", shotname, "rad.dat");
    int prof_d_length = read_vector(&prof_d, "input/ionProf", shotname, "ionyeald.dat");    
    
    int profx_r_length = read_vector(&profx_r, "input/ionProf", shotname, "xrad.dat", false);
    int profx_d_length = read_vector(&profx_d, "input/ionProf", shotname, "xionyeald.dat", false);
    
    if (prof_r_length != prof_d_length){
        printf("ERROR: Length of PROF_R and PROF_D are different!\n");
    }

    if (profx_r_length != profx_d_length){
        printf("ERROR: Length of PROFX_R and PROFX_D are different!\n");
    }

    if (profx_r_length <= 0){
        printf("Cross section beam profile: OFF\n");
    }else{
        printf("Cross section beam profile: ON\n");
    }

    double diam = beam.diameter / 1000.0;
    double deflH = beam.toroidal_deflection/180*PI;
    double deflV = beam.vertical_deflection/180*PI;
    
    Vabs = sqrt(2 * beam.energy*1000*ELEMENTARY_CHARGE/ beam.mass/ AMU);
    
    /* one-ion beam */
    if (profx_r_length > 0){
        for (i=0; i<profx_r_length; i++){
            profx_d[i] /= profx_r[i];
        }
    }
    /* initialize random generator */
    srand ( time(NULL) );
    for (i=0; i<shot.particle_number; i++){
        /* set position of particles */
        do{
            ionisation_yeald = (double)rand()/RAND_MAX;
            XR[i] = linear_interpolate(prof_d, prof_d_length, prof_r, prof_r_length, ionisation_yeald);
        }while (isnan(XR[i])||XR[i]<0);
        do{
            if (profx_r_length <= 0){
                XZ[i]=(double)(rand()-RAND_MAX/2)/RAND_MAX*diam;
                XT[i]=(double)(rand()-RAND_MAX/2)/RAND_MAX*diam;
            }else{
                ionisation_yeald = (double)rand()/RAND_MAX;
                xsec_ang = (double)rand()/RAND_MAX*2*PI;
                xsec_rad = linear_interpolate(profx_d, profx_d_length, profx_r, profx_r_length, ionisation_yeald)*(diam/2);
                XZ[i]= sin(xsec_ang) * xsec_rad;
                XT[i]= cos(xsec_ang) * xsec_rad;
            }
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
    if (x_length != y_length)   printf("ERROR: in interpolation. Two input vectors have different length.");    
    for (i=1; (i<x_length) && (x_vector[i-1]>x_value); i++);    
    if(i>1){--i;}else{i=1;}    
    return y_vector[i] - (y_vector[i]-y_vector[i-1])*(x_value-x_vector[i-1])/(x_vector[i]-x_vector[i-1]);
}

int get_array_size(double *array){
    return (int)(sizeof(array)/sizeof(array[0]));
} 
