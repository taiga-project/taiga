// TAIGA default parameters

#define $R_defl 2.3                 //! radial position of deflection plates in meter -> TOROIDAL DEFLECTION


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
    printf("%s\n", concat("TAIGA ", TAIGA_VERSION," (r", GIT_REV, ")"));
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
        CopyFile(run.parameter_file, concat(run.folder_out,"/",run.runnumber,"/parameters.sh"));
        
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
        init_coords(beam, shot, run, host_global, shared_global);
        
        //! grid
        init_grid(shot, run, host_common, shared_common);
        magnetic_field_read_and_init(shot, run, host_common, shared_common);
        if (shot.electric_field_on) shot.electric_field_on = electric_field_read_and_init(shot, run, host_common, shared_common);
        
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
        
        if (run.debug == 1 && !FASTMODE)   debug_message_init(host_global->rad, host_global->z, host_global->tor, host_global->vrad, host_global->vz, host_global->vtor);
        
        size_t dimX = host_global->particle_number*sizeof(double);
        
        init_device_structs(beam, shot, run, shared_global, shared_common);
        sync_device_structs(device_global, shared_global, device_common, shared_common);
        
        for (int step_i=0; step_i<run.step_host; ++step_i){
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
            
            if (run.debug == 1)    printf("Step\t%d/%d\n",step_i,run.step_host);
            // UNSOLVED: if (run.debug == 1 && !FASTMODE)    debug_message_run(host_global->rad, host_global->z, host_global->tor, host_global->vrad, host_global->vz, host_global->vtor);
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
        printf ("CPU->HDD copy time:  %lf s\n", run.cpu_time_copy);
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
        
        if (run.debug == 1)    debug_service_vars(host_service_array);
        
        //! CUDA profiler STOP
        cudaProfilerStop();
        
        fill_header_file(host_common, beam, shot, run);
        
        if (!FASTMODE){
            save_endpoints(host_global, run);
        }
        
        printf("\nData folder: %s/%s\n\n", run.folder_out, run.runnumber);
        
        //! Free CUDA
        /*for (int i=0; i<6; ++i){
            cudaFree(host_global.coords[i]);
        }
        cudaFree(host_global.coords);
        cudaFree(device_common.spline_grid[0]);
        cudaFree(device_common.spline_grid[1]);
        cudaFree(device_common.spline_grid);

        for (int i=0; i<3; ++i){
            cudaFree(host_common.espline[i]);
            cudaFree(host_common.bspline[i]);
        }
        cudaFree(host_common.espline);
        cudaFree(host_common.bspline);
        

        //! Free RAM    
        free(device_common.spline_grid[0]);
        free(device_common.spline_grid[1]);
        free(device_common.spline_grid);

        for (int i=0; i<3; ++i){
            free(host_common.espline[i]);
            free(host_common.bspline[i]);
        }
        free(host_common.espline);
        free(host_common.bspline);
        if (!FASTMODE){
            for (int i=0; i<6; ++i){
                free(host_global.coords[i]);
            }
            free(host_global.coords);
        }*/
        
        
        //! FREE host_service_array variables (RAM, cuda)
        free(host_service_array);  cudaFree(device_service_array);
        printf("Ready.\n\n");
    }
}

void print_run_details(TaigaGlobals *host_global, TaigaCommons *host_common, ShotProp shot, RunProp run){
    printf("===============================\n");
    printf("%s\n", concat("TAIGA ", TAIGA_VERSION," (r", GIT_REV, ")"));
    printf("Shotname: %s\n", shot.name); 
    printf("Detector: %s\n", shot.detector_mask);
    printf("  R:\t%lf\n", host_common->detector_geometry[0]);
    printf("  Z:\t%lf\n", host_common->detector_geometry[1]);
    printf("  T:\t%lf\n", host_common->detector_geometry[2]);
    printf("  angle (Z/R):\t%lf°\n", atan(host_common->detector_geometry[3])/PI*180.0);
    printf("  angle (T/R):\t%lf°\n", atan(host_common->detector_geometry[4])/PI*180.0);
    printf("===============================\n");
    printf("Number of blocks (threads): %d\n", run.block_number);
    printf("Block size: %d\n", run.block_size);
    printf("Number of particles: %d\n", host_global->particle_number);
    printf("Max steps on device (GPU): %d\n", run.step_device);
    printf("Max steps on host (HDD): %d\n", run.step_host);
    printf("===============================\n");
}

void fill_header_file(TaigaCommons *common, BeamProp beam, ShotProp shot, RunProp run){
    export_header(concat("TAIGA ", TAIGA_VERSION," (r", GIT_REV, ")"), run.folder_out, run.runnumber);
    export_header_addline(run.folder_out, run.runnumber);
    export_header(concat("Shot ID: ",shot.name), run.folder_out, run.runnumber);
    export_header(concat("Run ID:  ",run.runnumber), run.folder_out, run.runnumber);
    export_header_addline(run.folder_out, run.runnumber);
    export_header("ABP ION TRAJECTORIES", run.folder_out, run.runnumber);
    
    if(READINPUTPROF==1){
        export_header("Manual (6D) input profile", run.folder_out, run.runnumber);
    }else if(RENATE==110){
        export_header("TS + Renate 1.1.0 input profile", run.folder_out, run.runnumber);
    }
    export_header_addline(run.folder_out, run.runnumber);
    
    if(!READINPUTPROF){
        export_header("Beam energy", "keV", beam.energy, run.folder_out, run.runnumber);
        export_header("Atomic mass", "AMU", beam.mass, run.folder_out, run.runnumber);
        export_header("Beam diameter", "mm", beam.diameter*1000, run.folder_out, run.runnumber);
        export_header("Beam deflection (toroidal/vertical)", "°", beam.toroidal_deflection*180.0/PI, beam.vertical_deflection*180.0/PI, run.folder_out, run.runnumber);
    }
    
    export_header("Number of ions", "", (double)run.particle_number, run.folder_out, run.runnumber);
    export_header_addline(run.folder_out, run.runnumber);
    export_header("Detector position (R)", "m", common->detector_geometry[0], run.folder_out, run.runnumber);
    export_header("Detector position (Z)", "m", common->detector_geometry[1], run.folder_out, run.runnumber);
    export_header("Detector position (T)", "m", common->detector_geometry[2], run.folder_out, run.runnumber);
    export_header("Detector angle (Z/R)", "°", atan(common->detector_geometry[3])*180.0/PI, run.folder_out, run.runnumber);
    export_header("Detector angle (T/R)", "°", atan(common->detector_geometry[4])*180.0/PI, run.folder_out, run.runnumber);
    export_header(concat("Detector mask:  \t", shot.detector_mask), run.folder_out, run.runnumber);
    export_header_addline(run.folder_out, run.runnumber);
    export_header("Timestep", "s", run.timestep, run.folder_out, run.runnumber);
    export_header_addline(run.folder_out, run.runnumber);
    export_header("CUDA kernel runtime", "s", run.cuda_time_core, run.folder_out, run.runnumber);
    export_header("CUDA memcopy time", "s", run.cuda_time_copy, run.folder_out, run.runnumber);
    export_header("CPU->HDD copy time", "s", run.cpu_time_copy, run.folder_out, run.runnumber);
    export_header_addline(run.folder_out, run.runnumber);
    export_header("Number of blocks (threads)", "", run.block_number, run.folder_out, run.runnumber);
    export_header("Block size", "", run.block_size, run.folder_out, run.runnumber);
    export_header("Length of a loop", "", run.step_device, run.folder_out, run.runnumber);
    export_header("Number of loops", "", run.step_host, run.folder_out, run.runnumber);
}

void save_trajectories(TaigaGlobals *host_global, RunProp run){
    export_data(host_global->rad,  host_global->particle_number, run.folder_out, run.runnumber, "t_rad.dat");
    export_data(host_global->z,    host_global->particle_number, run.folder_out, run.runnumber, "t_z.dat");
    export_data(host_global->tor,  host_global->particle_number, run.folder_out, run.runnumber, "t_tor.dat");
    export_data(host_global->vrad, host_global->particle_number, run.folder_out, run.runnumber, "t_vrad.dat");
    export_data(host_global->vz,   host_global->particle_number, run.folder_out, run.runnumber, "t_vz.dat");
    export_data(host_global->vtor, host_global->particle_number, run.folder_out, run.runnumber, "t_vtor.dat");
}

void save_endpoints(TaigaGlobals *host_global, RunProp run){
    export_data(host_global->rad,  host_global->particle_number, run.folder_out, run.runnumber, "rad.dat");
    export_data(host_global->z,    host_global->particle_number, run.folder_out, run.runnumber, "z.dat");
    export_data(host_global->tor,  host_global->particle_number, run.folder_out, run.runnumber, "tor.dat");
    export_data(host_global->vrad, host_global->particle_number, run.folder_out, run.runnumber, "vrad.dat");
    export_data(host_global->vz,   host_global->particle_number, run.folder_out, run.runnumber, "vz.dat");
    export_data(host_global->vtor, host_global->particle_number, run.folder_out, run.runnumber, "vtor.dat");
    
    export_table(run.folder_out, run.runnumber, "coords.dat", host_global->particle_number,
        host_global->rad, "R [m]",      host_global->z, "Z [m]",      host_global->tor, "T [m]", 
        host_global->vrad, "v_R [m/s]", host_global->vz, "v_Z [m/s]", host_global->vtor, "v_T [m/s]");
}
