__global__ void detector_postproc(TaigaGlobals *global, TaigaCommons *common, DetectorProp *detector){
//(double **x_ptr, double *detector->xgrid, int detector->length_xgrid, double *detector->ygrid, int detector->length_ygrid, double *det, int *detcellid){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (global->detcellid[idx] == 0){
        double x, y;
        
        x = global->tor[idx] - common->detector_geometry[2];
        x*=1000; //in mm
        y = (common->detector_geometry[1] - global->z[idx]) / common->detector_geometry[3];
        y*=1000; //in mm
        
        int x_cellid = -1, y_cellid = -1;
        
        if ((x >= detector->xgrid[0]) & (x <= detector->xgrid[detector->length_xgrid*2-1]) & (y >= detector->ygrid[0]) & (y <= detector->ygrid[detector->length_ygrid*2-1])) {
            for (int i=0; i<(detector->length_xgrid); ++i) {
                if ((x >= detector->xgrid[2*i]) & (x <= detector->xgrid[2*i+1]))    x_cellid = i;
            }
            for (int j=0; j<detector->length_ygrid; ++j) {
                if ((y >= detector->ygrid[2*j]) & (y <= detector->ygrid[2*j+1]))    y_cellid = j;
            }
        }
        
        if ((x_cellid >= 0) & (y_cellid >= 0)) {
            global->detcellid[idx] = x_cellid * detector->length_ygrid + detector->length_ygrid-1 - y_cellid + 1;
        }
        global->detcellid[idx] =detector->xgrid[1];
    }
}
