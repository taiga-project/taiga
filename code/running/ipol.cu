__device__ double interpolate(double x1, double x2, double x, double y1, double y2, double x2_x1){ 
    return (x2-x)/x2_x1*y1 + (x-x1)/x2_x1*y2;
}

__device__ double interpolate(double x1, double x2, double x, double y1, double y2){
    return interpolate(x1, x2, y1, y2, x2-x1);
}

__device__ int calculate_detection_position(double position_rad, double position_z, double position_tor, double position_rad_prev, double position_z_prev, double position_tor_prev, double *detector_geometry, double *position, double speed_rad){

    int finished = 0;
    
    double detector_R = detector_geometry[0];
    double detector_z = detector_geometry[1];
    double detector_tan = detector_geometry[3];
    
    double detector_distance, detector_distance_prev;
    detector_distance = (position_rad-detector_R) + detector_tan*(position_z-detector_z);
    detector_distance_prev = (position_rad_prev-detector_R) + detector_tan*(position_z_prev-detector_z);
    
    if((detector_distance*detector_distance_prev<=0) && (speed_rad>0)){

        double detector_distance_change = (detector_distance-detector_distance_prev);
        position[0] = interpolate(position_rad_prev, position_rad, detector_distance_prev, detector_distance, detector_distance_change);
        position[1] = interpolate(position_z_prev,   position_z,   detector_distance_prev, detector_distance, detector_distance_change);
        position[2] = interpolate(position_tor_prev, position_tor, detector_distance_prev, detector_distance, detector_distance_change);
        finished = 1;
    }
    
    return finished;
}
