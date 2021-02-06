#define SERVICE_VAR_LENGTH 10

#define ERRORCHECK() cErrorCheck(__FILE__, __LINE__)

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <time.h>
#include <math.h>
#include <string.h>
//#include <filesystem>

#include <cuda_profiler_api.h>
//#include "test/cuda/nvToolsExt.h"

#include "taiga_constants.h"
#include "prop.h"
#include "main.cuh"
#include "taiga_write.cu"
#include "debug_functions.c"
#include "basic_functions.h"

#include "dataio/data_import.c"
#include "dataio/field_import.c"
#include "dataio/parameter_reader.c"

#include "taiga_init.c"
#include "dataio/beam.h"
#if READINPUTPROF == 1
    #include "dataio/beam_manual_profile.c"
#elif RENATE == 110
    #include "dataio/beam_renate110.c"
#else
    #error A valid beam module is required!
#endif

#include "dataio/data_export.c"

#include "running/rk4.cu"
#include "running/detection.cu"
#include "running/undetected.cu"
#include "running/cyl2tor.cu"
#include "running/traj.cu"
#include "running/generate_coords.cu"
#include "running/taiga.cu"

#include "detector_module.cu"
#include "running/detector_postproc.cu"
#include "running/detector_sum.cu"

void input_init_taiga(int argc, char *argv[], ShotProp *shot, BeamProp *beam, RunProp *run){
    
    char *input;
    for (int i=1; i<argc; ++i){
        input = strtok(argv[i], "=");
        if (!strcmp(input, "--debug") || !strcmp(input, "-d")){
            run->debug = 1;
        }else if (!strcmp(input, "--fulltrace") || !strcmp(input, "-f")){
            run->step_host = 2000;
            run->step_device = 1;
        }else if (!strcmp(input, "--help") || !strcmp(input, "-h")){
            run->help = 1;
        }else if (!strcmp(input, "--devices") || !strcmp(input, "-l")){   
            run->help = 2;
        }else if (!strcmp(input, "--parameter_file") || !strcmp(input, "-p")){
            input = strtok(NULL, "=");
            strcpy(run->parameter_file, input);
            printf("Parameter file: %s\n", run->parameter_file);
        }else if (!strcmp(input, "--runnumber_file")){
            strcpy(run->runnumber_file, input);
            printf("Runnumber file: %s\n", run->runnumber_file);
        }else if (!strcmp(input, "--runnumber") || !strcmp(input, "-r")){
            input = strtok(NULL, "=");
            strcpy(run->runnumber, input);
            strcpy(run->runnumber_file, "console init");
            printf("Runnumber: %s\n", run->runnumber);
        }else if (!strcmp(input, "--ion-source") || !strcmp(input, "-s")){
            input = strtok(NULL, "=");
            strcpy(run->ion_source_file, input);
            printf("Ion source file: %s\n", run->ion_source_file);
        }else if (!strcmp(input, "--ion-source-coords")){
            input = strtok(NULL, "=");
            strcpy(run->io_coordinate_order, input);
            printf("Order of coordinates in input file: %s\n", run->io_coordinate_order);
        }
    }
}

void print_help_message(){        
    printf("%s\n", concat("TAIGA ", TAIGA_VERSION," (r", GIT_REV, ")", NULL));
    printf("Usage: taiga.exe [options]\nOptions:\n");
    printf("  -d, --debug                 Print additional debug informations\n");
    printf("  -f, --fulltrace             Save coordinates at every timestep\n");
    printf("  -h, --help                  Help message\n");
    printf("  -l, --devices               List GPU devices\n");
    printf("  -p, --parameter_file=PATH   Parameter file path\n");
    printf("      --runnumber_file=PATH   Runnumber file path\n");
    printf("  -r  --runnumber=INTEGER     Runnumber value\n");
    printf("  -s, --ion-source=PATH       Ion source path\n");
    printf("      --ion-source-coords=XXX Order of coordinates (RZT or RTZ) in input file\n");
}

int main(int argc, char *argv[]){
    ShotProp shot; init_shot_prop(&shot);
    BeamProp beam; init_beam_prop(&beam);
    RunProp run;   init_run_prop(&run);
    input_init_taiga(argc, argv, &shot, &beam, &run);
    
    if (run.help == 1){
        print_help_message();
    }else if (run.help == 2){
        set_cuda(1);
    }else{
        parameter_reader(&beam, &shot, &run);
        runnumber_reader(&shot, &run);
        
        init_dir(run.folder_out, run.runnumber);
        CopyFile(run.parameter_file, concat(run.folder_out,"/",run.runnumber,"/parameters.sh", NULL));
        
        //! CUDA profiler START
        cudaProfilerStart();
        set_cuda(run.debug);
        
        TaigaGlobals *device_global, *host_global, *shared_global;
        TaigaCommons *device_common, *host_common, *shared_common;
        DetectorProp *shared_detector, *device_detector;
        
        size_t size_global = sizeof(TaigaGlobals);
        size_t size_commons = sizeof(TaigaCommons);
        size_t size_detector_prop = sizeof(DetectorProp);
        
        host_global = (TaigaGlobals*)malloc(size_global);
        shared_global = (TaigaGlobals*)malloc(size_global);
        host_common = (TaigaCommons*)malloc(size_commons);
        shared_common = (TaigaCommons*)malloc(size_commons);
        shared_detector = (DetectorProp*)malloc(size_detector_prop);
        
        cudaMalloc((void **) &device_global, size_global);
        cudaMalloc((void **) &device_common, size_commons);
        cudaMalloc((void **) &device_detector, size_detector_prop);
        
        init_host(host_global, host_common);
        
        set_particle_number(&run, host_global, shared_global);
        
        //! coordinates
        init_coords(&beam, &shot, &run, host_global, shared_global);
        
        //! grid
        init_grid(shot, run, host_common, shared_common);
        magnetic_field_read_and_init(shot, run, host_common, shared_common);
        if (shot.is_electric_field_on) shot.is_electric_field_on = electric_field_read_and_init(shot, run, host_common, shared_common);
        
        // detector
        set_detector_geometry(shot, host_common, shared_common);
        init_detector(shared_detector, device_detector, shot);
        
        // <service value>
        size_t dimService = SERVICE_VAR_LENGTH * sizeof(double);
        double *host_service_array, *device_service_array;
        host_service_array = (double *)malloc(dimService);
        
        for(int i=0; i<SERVICE_VAR_LENGTH; ++i){
            host_service_array[i] = 0;
        }
        
        host_service_array[4] = 55555.55555;
        cudaMalloc((void **) &device_service_array,  dimService);
        cudaMemcpy(device_service_array, host_service_array, dimService, cudaMemcpyHostToDevice);
        // </service value>
        
        if (!FASTMODE){
           save_trajectories(host_global, run);
        }
        
        print_run_details(host_global, host_common, shot, run);
        
        //! Set CUDA timer 
        cudaEvent_t cuda_event_core_start, cuda_event_core_end, cuda_event_copy_start, cuda_event_copy_end;
        clock_t cpu_event_copy_start, cpu_event_copy_end;
        float cuda_event_core, cuda_event_copy;
        cudaEventCreate(&cuda_event_core_start);
        cudaEventCreate(&cuda_event_core_end);
        cudaEventCreate(&cuda_event_copy_start);
        cudaEventCreate(&cuda_event_copy_end);
        
        if (run.debug == 1 && !FASTMODE)   debug_message_init(host_global);
        
        size_t dimX = host_global->particle_number*sizeof(double);
        
        init_device_structs(beam, shot, run, shared_global, shared_common);
        sync_device_structs(device_global, shared_global, device_common, shared_common);
        if (FASTMODE)   init_fastmode(beam, shot, run, device_global);
        
        for (long step_i=0; step_i<run.step_host; ++step_i){
            if (step_i == 0) cudaEventRecord(cuda_event_core_start, 0);
            
            taiga <<< run.block_number, run.block_size >>> (device_global, device_common, device_service_array);
            
            if (step_i == 0) cudaEventRecord(cuda_event_core_end, 0);
            cudaEventSynchronize(cuda_event_core_end);
            //ERRORCHECK();
            
            if (!FASTMODE){
                // ION COORDS (device2HOST)
                if (step_i == 0) cudaEventRecord(cuda_event_copy_start, 0);
                coord_memcopy_back(beam, shot, run, host_global, shared_global);
                //ERRORCHECK();
                if (step_i == 0) cudaEventRecord(cuda_event_copy_end, 0);
                
                // Save data to files
                cpu_event_copy_start = clock();
                save_trajectories(host_global, run);
                cpu_event_copy_end = clock();
            }
            
            if (run.debug == 1)    printf("Step\t%ld/%ld\n",step_i,run.step_host);
            if (run.debug == 1 && !FASTMODE)    debug_message_run(host_global);
        }
        
        // Get CUDA timer 
        cudaEventElapsedTime(&cuda_event_core, cuda_event_core_start, cuda_event_core_end);
        cudaEventElapsedTime(&cuda_event_copy, cuda_event_copy_start, cuda_event_copy_end);
        if (!FASTMODE) run.cpu_time_copy = ((double) (4.0+run.step_host)*(cpu_event_copy_end - cpu_event_copy_start)) / CLOCKS_PER_SEC;
        run.cuda_time_copy = (double) (1.0+run.step_host)*cuda_event_copy/1000.0;
        run.cuda_time_core =  run.step_host*cuda_event_core/1000.0;
        printf("===============================\n");
        printf ("CUDA kernel runtime: %lf s\n", run.cuda_time_core);
        printf ("CUDA memcopy time:   %lf s\n", run.cuda_time_copy);
        if (!FASTMODE)  printf ("CPU->HDD copy time:  %lf s\n", run.cpu_time_copy);
        printf("===============================\n");
        
        //UNSOLVED: undetected <<<1,1>>>(detcellid, host_global->particle_number, device_service_array);
        //UNSOLVED: printf("Lost particle ratio: \t %.4lf % \n\n", host_service_array[1]*100);
        
        //! MEMCOPY (device2HOST)
        cudaMemcpy(host_service_array, device_service_array, dimService, cudaMemcpyDeviceToHost);
        if(host_service_array[0] != 42.24){
            printf("\n +----------------------------+\n | Fatal error in running.    | \n | The CUDA did not run well. |\n | Service value: %11lf |\n +----------------------------+\n\n", host_service_array[0]);
        }else{
            printf("\nSuccessful run. \n\n");
        }
        
        detector_postproc <<< run.block_number, run.block_size >>> (device_global, device_common, device_detector);
        detector_sum <<<1,1>>> (device_global, device_common, device_detector);
        export_detector(shared_detector, device_detector, shared_global, shot, run);
        
        //! CUDA profiler STOP
        cudaProfilerStop();
        
        if (run.debug == 1)    debug_service_vars(host_service_array);
        
        fill_header_file(host_common, beam, shot, run);
        
        if (!FASTMODE){
            save_endpoints(host_global, run);
        }
        
        printf("\nData folder: %s/%s\n\n", run.folder_out, run.runnumber);
        
        //! FREE host_service_array variables (RAM, cuda)
        free(host_service_array);  cudaFree(device_service_array);
        printf("Ready.\n\n");
    }
}
