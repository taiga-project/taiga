__device__ double hermite_spline(double y1, double y2, double vy1, double vy2, double t){ 
    return (2*t*t*t-3*t*t+1)*x1 + (t*t*t-2*t*t+t)*v1 +(-2*t*t*t+3*t*t)*x2 + (t*t*t-t*t)*v2; 
}

__device__ double interpolate(double y1, double y2, double x, double x1, double x2, double x2_x1){ 
    return (x2-x)/x2_x1*y1 + (x-x1)/x2_x1*y2;
}

__device__ double interpolate(double y1, double y2, double x, double x1, double x2){
    return interpolate(y1, y2, x, x1, x2, x2-x1);
}

__device__ int calculate_detection_position(double *X, double *X_prev, double *detector_geometry){

    int finished = 0;
    
    double detector_R   = detector_geometry[0];
    double detector_z   = detector_geometry[1];
    double detector_tan = detector_geometry[3];
    
    double detector_distance = (X[0]-detector_R) + detector_tan*(X[1]-detector_z);
    double detector_distance_prev = (X_prev[0]-detector_R) + detector_tan*(X_prev[1]-detector_z);
    
    if((detector_distance*detector_distance_prev<=0) && (X[3]>0)){
        double detector_cos = 1/sqrt(1+detector_tan*detector_tan);
        detector_distance *= detector_cos;
        detector_distance_prev *= detector_cos;

        double X_new[6];     
        double detector_distance_rate = -detector_distance_prev/(detector_distance-detector_distance_prev);
        
        X_new[0] = hermite_spline(X_prev[0], X[0], X_prev[3], X[3], detector_distance_rate);
        X_new[1] = hermite_spline(X_prev[1], X[1], X_prev[4], X[4], detector_distance_rate);
        X_new[2] = hermite_spline(X_prev[2], X[2], X_prev[5], X[5], detector_distance_rate);
        X_new[3] = interpolate(X_prev[3], X[3], 0, detector_distance_prev, detector_distance);
        X_new[4] = interpolate(X_prev[4], X[4], 0, detector_distance_prev, detector_distance);
        X_new[5] = interpolate(X_prev[5], X[5], 0, detector_distance_prev, detector_distance);
        
        X[0] = X_new[0];
        X[1] = X_new[1];
        X[2] = X_new[2];
        X[3] = X_new[3];
        X[4] = X_new[4];
        X[5] = X_new[5];
        
        finished = 1;
        
    }
    
    return finished;
}
