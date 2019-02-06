__global__ void detector_postproc(double **x_ptr, double *det_x, int N_det_x, double *det_y, int N_det_y, double *det, int *detcellid){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (detcellid[idx] == 0){
        double x,y;
        
        x = x_ptr[2][idx] - det[2];                x*=1000; //in mm
        y = (det[1] - x_ptr[1][idx]) / det[3];     y*=1000; //in mm
        
        int x_cellid = -1, y_cellid = -1;
        
        if ((x >= det_x[0]) & (x <= det_x[N_det_x-1]) & (y >= det_y[0]) & (y <= det_y[N_det_y-1])) {
            for (int i=0; i<(N_det_x/2); i++) {
                if ((x >= det_x[2*i]) & (x <= det_x[2*i+1]))    x_cellid = i;
            }
            for (int j=0; j<N_det_y/2; j++) {
                if ((y >= det_y[2*j]) & (y <= det_y[2*j+1]))    y_cellid = j;
            }
        }
        
        if ((x_cellid >= 0) & (y_cellid >= 0)) {
            detcellid[idx] = (N_det_x/2-1 - x_cellid) * N_det_y/2 +  N_det_y/2-1 - y_cellid + 1;
        }
    }
}
