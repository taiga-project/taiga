#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>

void export_detector(DetectorProp* shared_detector, DetectorProp *device_detector, TaigaGlobals *shared_global, ShotProp shot, RunProp run){
    if (shared_detector->detector_module_on){
        int *host_counter;
        size_t size_counter = shared_detector->number_of_detector_cells * sizeof(int);
        
        host_counter = (int*)malloc(size_counter);
        cudaMemcpy(host_counter, shared_detector->counter, size_counter, cudaMemcpyDeviceToHost);
        export_data(host_counter, shared_detector->number_of_detector_cells, run.folder_out, run.runnumber, "detector", "cellcounter.dat", shared_detector->length_ygrid);
        
        if (!FASTMODE){
            size_t size_detcellid = run.block_size * run.block_number * sizeof(int);
            int* host_detcellid =(int*)malloc(size_detcellid);
            cudaMemcpy(host_detcellid, shared_global->detcellid, size_detcellid, cudaMemcpyDeviceToHost);
            export_data(host_detcellid, shared_global->particle_number, run.folder_out, run.runnumber, "detector", "cellid.dat");
        }
    }
    
    CopyFile(concat("input/detector/", shot.detector_mask, "/detx"), concat(run.folder_out,"/",run.runnumber,"/detector/detx"));
    
    CopyFile(concat("input/detector/", shot.detector_mask, "/dety"), concat(run.folder_out,"/",run.runnumber,"/detector/dety"));
}
