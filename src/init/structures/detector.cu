void init_detector(DetectorProp* shared_detector, DetectorProp *device_detector, ShotProp shot){
    double *host_detector_xgrid, *device_detector_xgrid;
    double *host_detector_ygrid, *device_detector_ygrid;
    double *device_counter;
    
    shared_detector->length_xgrid = read_vector(&host_detector_xgrid, "input/detector", shot.detector_mask, "detx", false)/2;
    shared_detector->length_ygrid = read_vector(&host_detector_ygrid, "input/detector", shot.detector_mask, "dety", false)/2;
    shared_detector->number_of_detector_cells = shared_detector->length_xgrid * shared_detector->length_ygrid;

    shared_detector->detector_module_on = (( shared_detector->length_xgrid > 0) && ( shared_detector->length_ygrid > 0));
    if (shared_detector->detector_module_on) {
        printf("===============================\n");
        printf("Detector postprocessor module: ON (%d x %d / %s)\n", shared_detector->length_xgrid, shared_detector->length_ygrid, shot.detector_mask);
        size_t size_detector_xgrid = 2 * shared_detector->length_xgrid * sizeof(double);
        size_t size_detector_ygrid = 2 * shared_detector->length_ygrid * sizeof(double);
        size_t size_counter = shared_detector->number_of_detector_cells * sizeof(double);
        size_t size_detector = sizeof(DetectorProp);
        
        cudaMalloc((void **) &device_detector_xgrid, size_detector_xgrid);
        cudaMalloc((void **) &device_detector_ygrid, size_detector_ygrid);
        cudaMalloc((void **) &device_counter, size_counter);
        
        cudaMemcpy(device_detector_xgrid, host_detector_xgrid, size_detector_xgrid, cudaMemcpyHostToDevice);
        cudaMemcpy(device_detector_ygrid, host_detector_ygrid, size_detector_ygrid, cudaMemcpyHostToDevice);
        shared_detector->xgrid = device_detector_xgrid;
        shared_detector->ygrid = device_detector_ygrid;
        shared_detector->counter = device_counter;
        
        cudaMemcpy(device_detector, shared_detector, size_detector, cudaMemcpyHostToDevice);
    }else{
        printf("===============================\n");
        printf("Detector postprocessor module: OFF\n");
    }
}

void set_detector_geometry(ShotProp shot, TaigaCommons *host_common, TaigaCommons *shared_common, DetectorProp *shared_detector){
    if ( (strcmp(shot.detector_mask, "no") == 0) ||
         (strcmp(shot.detector_mask, "no_mask") == 0) ||
         (strcmp(shot.detector_mask, "no mask") == 0)) {
        shared_detector->detector_module_on = false;
    } else {
        char *tokaniser;
        double *shared_detector_geometry;
        shared_detector->detector_module_on = true;
        size_t size_detector = 5 * sizeof(double);
        host_common->detector_geometry = (double *) malloc(size_detector);
        tokaniser = strtok(shot.detector_geometry, ",");
        host_common->detector_geometry[0] = strtod(tokaniser, NULL);
        tokaniser = strtok(NULL, ",");
        host_common->detector_geometry[1] = strtod(tokaniser, NULL);
        tokaniser = strtok(NULL, ",");
        host_common->detector_geometry[2] = strtod(tokaniser, NULL);
        tokaniser = strtok(NULL, ",");
        host_common->detector_geometry[3] = (strtod(tokaniser, NULL) * PI / 180.0);
        tokaniser = strtok(NULL, ",");
        host_common->detector_geometry[4] = (strtod(tokaniser, NULL) * PI / 180.0);

        cudaMalloc((void **) &shared_detector_geometry, size_detector);
        cudaMemcpy(shared_detector_geometry, host_common->detector_geometry, size_detector, cudaMemcpyHostToDevice);
        shared_common->detector_geometry = shared_detector_geometry;
    }
}