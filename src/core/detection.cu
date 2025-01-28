#include "core/maths/maths.cuh"

__device__ int calculate_detection_position(double *X, double *X_prev, double *detector_geometry, double timestep,
                                            int method){
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
        if (method == BEZIER_INTERPOLATION) {
            double detector_plane[4] = {cos_alpha, sin_alpha * cos_beta, sin_beta,
                                        -(cos_alpha * detector_R +
                                        sin_alpha * cos_beta * detector_Z +
                                        sin_beta * detector_T)};
            interpolate_bezier(X_prev, X, detector_plane, timestep);
        }else{
            for (int i = 0; i < 6; ++i)
                X[i] = interpolate(X_prev[i], X[i], 0, detector_distance_prev, detector_distance);
        }
        X[TIME_OF_FLIGHT_ID] += interpolate(0, 1, 0, detector_distance_prev, detector_distance) * timestep;
        return CALCULATION_FINISHED;
    }
    X[TIME_OF_FLIGHT_ID] += timestep;
    return CALCULATION_NOT_FINISHED;
}
