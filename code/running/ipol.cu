__device__ int calculate_detection_position(double position_rad, double position_z, double position_tor, double position_rad_prev, double position_z_prev, double position_tor_prev, double *detector_geometry, double *position, double speed_rad){

    int finished = 0;
    
    double detector_R = detector_geometry[0];
    double detector_z = detector_geometry[1];
    double detector_tan = detector_geometry[3];
    
    double detector_distance, detector_distance_prev;
    detector_distance = (position_rad-detector_R) + detector_tan*(position_z-detector_z);
    detector_distance_prev = (position_rad_prev-detector_R) + detector_tan*(position_z_prev-detector_z);
    
    if((detector_distance*detector_distance_prev<=0) && (speed_rad>0)){

        double di_rate = (-detector_distance_prev)/(detector_distance-detector_distance_prev);
        position[0] = position_rad_prev + di_rate*(position_rad-position_rad_prev);;
        position[1] = position_z_prev   + di_rate*(position_z-position_z_prev);
        position[2] = position_tor_prev + di_rate*(position_tor-position_tor_prev);
        finished = 1;
    }
    
    return finished;
}
