#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include "running/detector_postproc.cu"
#include "running/detector_sum.cu"

//void detector_module(double **x_ptr, double *detector, int *detcellid, char *shot.detector_mask, int max_blocks, int shot_block_size, int number_of_particles, char *export_folder, char *runnumber){

void init_detector(DetectorProp* shared_detector, DetectorProp *device_detector, ShotProp shot){
    double *device_detector_xgrid;
    double *device_detector_ygrid;
    int *device_counter;
    
    shared_detector->length_xgrid = read_vector(&shared_detector->xgrid, "input/detector", shot.detector_mask, "detx", false)/2;
    shared_detector->length_ygrid = read_vector(&shared_detector->ygrid, "input/detector", shot.detector_mask, "dety", false)/2;
    shared_detector->number_of_detector_cells = shared_detector->length_xgrid * shared_detector->length_ygrid;

    shared_detector->detector_module_on = (( shared_detector->length_xgrid > 0) && ( shared_detector->length_ygrid > 0));
    if (shared_detector->detector_module_on) {
        printf("===============================\n");
        printf("Detector postprocessor module: ON (%d x %d)\n", shared_detector->length_xgrid, shared_detector->length_ygrid);
        size_t size_detector_xgrid = 2 * shared_detector->length_xgrid * sizeof(double);
        size_t size_detector_ygrid = 2 * shared_detector->length_ygrid * sizeof(double);
        size_t size_counter = shared_detector->number_of_detector_cells * sizeof(int);
        
        cudaMalloc((void **) &device_detector_xgrid,  size_detector_xgrid);
        cudaMalloc((void **) &device_detector_ygrid,  size_detector_ygrid);
        cudaMalloc((void **) &device_counter,  size_counter);
        
        cudaMemcpy(device_detector_xgrid, shared_detector->xgrid, size_detector_xgrid, cudaMemcpyHostToDevice);
        cudaMemcpy(device_detector_ygrid, shared_detector->ygrid, size_detector_ygrid, cudaMemcpyHostToDevice);
        cudaMemcpy(device_counter, shared_detector->counter, size_counter, cudaMemcpyHostToDevice);
    }else{
        printf("===============================\n");
        printf("Detector postprocessor module: OFF\n");
    }
}

void export_detector(DetectorProp* shared_detector, DetectorProp *device_detector, TaigaGlobals *device_global, TaigaCommons *device_common, ShotProp shot, RunProp run){
    if (shared_detector->detector_module_on){
        int *host_counter;
        size_t size_counter = shared_detector->number_of_detector_cells * sizeof(int);
        
        host_counter = (int*)malloc(size_counter);
        cudaMemcpy(host_counter, shared_detector->counter, size_counter, cudaMemcpyDeviceToHost);
        export_data(host_counter, shared_detector->number_of_detector_cells, run.folder_out, run.runnumber, "detector", "cellcounter.dat", shared_detector->length_ygrid);
        
        /*printf("number of cells:%d\n",shared_detector->number_of_detector_cells);
        for(int i=0; i<shared_detector->number_of_detector_cells; ++i){
            printf("%2d. %d\n",i, shared_detector->counter[i]);
        }*/
    }
    
    CopyFile(concat("input/detector/", shot.detector_mask, "/detx"), concat(run.folder_out,"/",run.runnumber,"/detector/detx"));
    
    CopyFile(concat("input/detector/", shot.detector_mask, "/dety"), concat(run.folder_out,"/",run.runnumber,"/detector/dety"));
        
}
