__global__ generate_coords(int idx, double **position_all, double ***speed_all, double eperm, int *prof_size, double *prof_x, double *prof_d, double *profx_r, double *profx_d){
    int i;
    double Vabs, ionisation_yeald, xsec_rad, xsec_ang;
    
    Vabs = sqrt(2 * eperm);
    
    /* cross section normalisation */
    if (prof_size[1] > 0){
        for (i=0; i<prof_size[1]; i++){
            profx_d[i] /= profx_r[i];
        }
    }
    
    /* initialize random generator */
    srand ( time(NULL) );
    
    /* set position of particles */
    do{
        ionisation_yeald = (double)rand()/RAND_MAX;
        position_all[0][idx] = linear_interpolate(prof_d, prof_size[0], prof_r, prof_size[0], ionisation_yeald);
    }while (isnan(XR[i])||XR[i]<0);
    do{
        //if (prof_size[1] <= 0){
            position_all[1][idx]=(double)(rand()-RAND_MAX/2)/RAND_MAX*beam.diameter;
            position_all[2][idx]=(double)(rand()-RAND_MAX/2)/RAND_MAX*beam.diameter;
        /*}else{
            ionisation_yeald = (double)rand()/RAND_MAX;
            xsec_ang = (double)rand()/RAND_MAX*2*PI;
            xsec_rad = linear_interpolate(profx_d, prof_size[1], profx_r, prof_size[1], ionisation_yeald)*(beam.diameter/2);
            XZ[i]= sin(xsec_ang) * xsec_rad;
            XT[i]= cos(xsec_ang) * xsec_rad;
        }*/
    }while ((XZ[i]*XZ[i]+XT[i]*XT[i])>=(beam.diameter/2)*(beam.diameter/2));
    
    /* toroidal deflection */
    position_all[2][idx] += tan(beam.toroidal_deflection) * ($R_defl - XR[i]);
    
    /* set velocity of particles */
    speed_all[0][idx] = -Vabs*cos(beam.vertical_deflection)*cos(beam.toroidal_deflection);
    speed_all[1][idx] =  Vabs*sin(beam.vertical_deflection);
    speed_all[2][idx] =  Vabs*cos(beam.vertical_deflection)*sin(beam.toroidal_deflection);
    
}
