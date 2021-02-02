__device__ __forceinline__ double calculate_detector_x(double particle_tor, double detector_tor){
    return (particle_tor - detector_tor)*1000.0; //in mm
}

__device__ __forceinline__ double calculate_detector_y(double particle_z, double detector_z, double angle){
    return (particle_z - detector_z)/cos(angle)*1000.0; //in mm
}

__device__ int get_cell_array_index(double value, double *array, int cell_number){
    for (int i=0; i<(cell_number); ++i) {
        if ((value >= array[2*i]) & (value <= array[2*i+1]))    return i;
    }
    return OUT_OF_RANGE;
}

__device__ __forceinline__ bool is_in_range(double value, double limit1, double limit2){
    return ( ( (value >=limit1)&(value<=limit2) ) | ( (value >=limit2)&(value<=limit1) ) );
}

__global__ void detector_postproc(TaigaGlobals *global, TaigaCommons *common, DetectorProp *detector){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (global->detcellid[idx] == CALCULATION_FINISHED){
        double x, y;
        x = calculate_detector_x(global->tor[idx], common->detector_geometry[2]);
        y = calculate_detector_y(global->z[idx], common->detector_geometry[1], common->detector_geometry[3]);
        
        int x_cellid = OUT_OF_RANGE, y_cellid = OUT_OF_RANGE;
        
        if ( is_in_range(x, detector->xgrid[0], detector->xgrid[detector->length_xgrid*2-1]) & 
             is_in_range(y, detector->ygrid[0], detector->ygrid[detector->length_ygrid*2-1])){
            x_cellid = get_cell_array_index(x, detector->xgrid, detector->length_xgrid);
            y_cellid = get_cell_array_index(y, detector->ygrid, detector->length_ygrid);
        }
        
        if ((x_cellid != OUT_OF_RANGE) & (y_cellid != OUT_OF_RANGE)) {
            global->detcellid[idx] = x_cellid * detector->length_ygrid + y_cellid + 1;
        }
    }
}
