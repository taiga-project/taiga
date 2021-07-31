#include "maths.cuh"

__device__ int calculate_detection_position(double *X, double *X_prev, double *detector_geometry){
    double detector_R   = detector_geometry[0];
    double detector_Z   = detector_geometry[1];
    double detector_T   = detector_geometry[2];
    double sin_alpha = sin(detector_geometry[3]);
    double cos_alpha = cos(detector_geometry[3]);
    double sin_beta = sin(detector_geometry[4]);
    double cos_beta = cos(detector_geometry[4]);
    
    double detector_distance =
        ( (X[0]-detector_R)*cos_alpha +
        (X[1]-detector_Z)*sin_alpha )*cos_beta + 
        (X[2]-detector_T)*sin_beta;
    double detector_distance_prev = 
        ( (X_prev[0]-detector_R)*cos_alpha +
        (X_prev[1]-detector_Z)*sin_alpha )*cos_beta +
        (X_prev[2]-detector_T)*sin_beta;
    
    if((detector_distance*detector_distance_prev<=0) && (X[3]>0)){
        double X_new[6];
        double detector_distance_rate = -detector_distance_prev/(detector_distance-detector_distance_prev);
        
        for (int i=0; i<6; ++i) X[i] = interpolate(X_prev[i], X[i], 0, detector_distance_prev, detector_distance);
        
        return CALCULATION_FINISHED;
    }
    
    return CALCULATION_NOT_FINISHED;
}
