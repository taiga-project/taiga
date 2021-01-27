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

#include "detector_module.c"
//#include "running/detector_postproc.cu"

void input_init_taiga(int argc, char *argv[], ShotProp *shot, BeamProp *beam, RunProp *run){
    
    char *input;
    for (int i=1; i<argc; i++){
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
        TaigaGlobals *dev_global, *host_global, *shared_global;
        TaigaCommons *dev_common, *host_common, *shared_common;
        
        size_t dim_global = sizeof(TaigaGlobals);
        size_t dim_commons = sizeof(TaigaCommons);
        
        host_global = (TaigaGlobals*)malloc(dim_global);
        shared_global = (TaigaGlobals*)malloc(dim_global);
        host_common = (TaigaCommons*)malloc(dim_commons);
        shared_common = (TaigaCommons*)malloc(dim_commons);
        
        cudaMalloc((void **) &dev_global, dim_global);
        cudaMalloc((void **) &dev_common, dim_commons);
        
        init_host(host_global, host_common);
        
        parameter_reader(&shot, &beam, &run);
        runnumber_reader(&shot, &run);
        
        init_dir(run.folder_out, run.runnumber);
        CopyFile(run.parameter_file, concat(run.folder_out,"/",run.runnumber,"/parameters.sh"));

        size_t dimD = 5 * sizeof(double);
        double *DETECTOR, *detector;
        DETECTOR = (double *)malloc(dimD);  cudaMalloc((void **) &detector,  dimD);
        
        set_detector_geometry(DETECTOR, shot.detector_geometry);
        
        printf("%s\n", concat("TAIGA ", TAIGA_VERSION," (r", GIT_REV, ")"));
        printf("Shotname: %s\n", shot.name); 
        printf("Detector: %s\n", shot.detector_mask);
        printf("  R:\t%lf\n", DETECTOR[0]);
        printf("  Z:\t%lf\n", DETECTOR[1]);
        printf("  T:\t%lf\n", DETECTOR[2]);
        printf("  angle (Z/R):\t%lf°\n", atan(DETECTOR[3])/PI*180.0);
        printf("  angle (T/R):\t%lf°\n", atan(DETECTOR[4])/PI*180.0);
        printf("===============================\n");
        
        host_global->particle_number = run.block_size * run.block_number;
        
        if (READINPUTPROF == 1){
            double *X_temp;
            host_global->particle_number = read_vector(&X_temp, "input", "manual_profile", "rad.dat");
            run.block_number = host_global->particle_number / run.block_size+1;
            free(X_temp);
        }
        
        //! CUDA profiler START
        cudaProfilerStart();
        
        set_cuda(run.debug);
        
        // set timestamp
        time_t rawtime;
        struct tm *info;
        
        printf("Number of blocks (threads): %d\n", run.block_number);
        printf("Block size: %d\n", run.block_size);
        printf("Number of particles: %d\n", host_global->particle_number);
        printf("Max steps on device (GPU): %d\n", run.step_device);
        printf("Max steps on host (HDD): %d\n", run.step_host);
        
        //! coordinates
        init_coords(beam, shot, run, host_global, shared_global);
        
        //! grid
        init_grid(shot, run, host_common, shared_common);
        int magnetic_field_loaded = magnetic_field_read_and_init(shot, run, host_common, shared_common);
        if (shot.electric_field_module) shot.electric_field_module = electric_field_read_and_init(shot, run, host_common, shared_common);
        
        // detector cell id
        size_t size_detcellid = host_global->particle_number * sizeof(int);
        int *DETCELLID, *detcellid;
        DETCELLID = (int *)malloc(size_detcellid); cudaMalloc((void **) &detcellid, size_detcellid);
        
        // service value
        size_t dimService = SERVICE_VAR_LENGTH * sizeof(double);
        double *SERVICE_VAR, *service_var;
        SERVICE_VAR = (double *)malloc(dimService);
        
        for(int i=0 ; i<SERVICE_VAR_LENGTH ; ++i){
            SERVICE_VAR[i] = 0;
        }
        
        SERVICE_VAR[4] = 55555.55555;
        cudaMalloc((void **) &service_var,  dimService);
        cudaMemcpy(service_var, SERVICE_VAR, dimService, cudaMemcpyHostToDevice);
        printf("temp.\n");
        
        //! MEMCOPY (HOST2device)
        printf("L209\n");
        //! DETECTOR COORDS (HOST2device)
        shared_common->detector_geometry = DETECTOR;
        printf("L212\n");
        size_t size_commons = sizeof(TaigaCommons);
        
        //cudaMemcpy(dev_common, shared_common, size_commons, cudaMemcpyHostToDevice);
        
        printf("L217\n");
        if (!FASTMODE){
        printf("L219\n");
            // OUTPUT INIT
            export_data(host_global->rad,  host_global->particle_number, run.folder_out, run.runnumber, "t_rad.dat");
            export_data(host_global->z,    host_global->particle_number, run.folder_out, run.runnumber, "t_z.dat");
            export_data(host_global->tor,  host_global->particle_number, run.folder_out, run.runnumber, "t_tor.dat");
            export_data(host_global->vrad, host_global->particle_number, run.folder_out, run.runnumber, "t_vrad.dat");
            export_data(host_global->vz,   host_global->particle_number, run.folder_out, run.runnumber, "t_vz.dat");
            export_data(host_global->vtor, host_global->particle_number, run.folder_out, run.runnumber, "t_vtor.dat");
        }
        
        printf("L229\n");
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
        
        //dev_common.step_counter = 0;
        printf("L244\n");
        init_device_structs(beam, shot, run, dev_global, shared_global, dev_common, shared_common);
        printf("L247\n");
        
        for (int step_i=0;step_i<run.step_host;step_i++){
            
            if (step_i == 0) cudaEventRecord(cuda_event_core_start, 0);
           printf("L250\n");
            taiga <<< run.block_number, run.block_size >>> (dev_global, dev_common, service_var);
            printf("L252\n");
            if (step_i == 0) cudaEventRecord(cuda_event_core_end, 0);
            cudaEventSynchronize(cuda_event_core_end);
            //ERRORCHECK();
            
            if (!FASTMODE){
                // ION COORDS (device2HOST)
                if (step_i == 0) cudaEventRecord(cuda_event_copy_start, 0);
                //coord_memcopy_back(host_global, shared_global, dev_common, beam, shot, run);
                //ERRORCHECK();                
                if (step_i == 0) cudaEventRecord(cuda_event_copy_end, 0);
                
                // Save data to files
                cpu_event_copy_start = clock();/*
                export_data(host_global->rad,  host_global->particle_number, run.folder_out, run.runnumber, "t_rad.dat");
                export_data(host_global->z,    host_global->particle_number, run.folder_out, run.runnumber, "t_z.dat");
                export_data(host_global->tor,  host_global->particle_number, run.folder_out, run.runnumber, "t_tor.dat");
                export_data(host_global->vrad, host_global->particle_number, run.folder_out, run.runnumber, "t_vrad.dat");
                export_data(host_global->vz,   host_global->particle_number, run.folder_out, run.runnumber, "t_vz.dat");
                export_data(host_global->vtor, host_global->particle_number, run.folder_out, run.runnumber, "t_vtor.dat");*/
                cpu_event_copy_end = clock();
            }
            
            if (run.debug == 1)    printf("Step\t%d/%d\n",step_i,run.step_host);
            if (run.debug == 1 && !FASTMODE)    debug_message_run(host_global->rad, host_global->z, host_global->tor, host_global->vrad, host_global->vz, host_global->vtor);
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
        
        undetected <<<1,1>>>(detcellid, host_global->particle_number, service_var);
        
        //! MEMCOPY (device2HOST)
        
        cudaMemcpy(SERVICE_VAR, service_var, dimService, cudaMemcpyDeviceToHost);
        if(SERVICE_VAR[0]!=42.24){
            printf("\n +----------------------------+\n | Fatal error in running.    | \n | The CUDA did not run well. |\n | Service value: %11lf |\n +----------------------------+\n\n", SERVICE_VAR[0]);
        }else{
            printf("\nSuccessful run. \n\n");
        }
        
        printf("Lost particle ratio: \t %.4lf % \n\n", SERVICE_VAR[1]*100);
        
        /* --> detector_module(x_ptr, detector, detcellid, shot.detector_mask, run.block_number, run.block_size, host_global->particle_number, run.folder_out, run.runnumber);
        cudaMemcpy(DETCELLID, detcellid, size_detcellid, cudaMemcpyDeviceToHost);
        export_data(DETCELLID, host_global->particle_number, run.folder_out, run.runnumber, "detector", "cellid.dat"); <--*/
        
        if (run.debug == 1)    debug_service_vars(SERVICE_VAR);
        
        //! CUDA profiler STOP
        cudaProfilerStop();
        
        fill_header_file(shot, beam, run, DETECTOR);
        
        /*if (!FASTMODE){
            //! Save data to files
            export_data(host_global.rad, host_global.particle_number, run.folder_out, run.runnumber, "rad.dat");
            export_data(host_global.z,   host_global.particle_number, run.folder_out, run.runnumber, "z.dat");
            export_data(host_global.tor, host_global.particle_number, run.folder_out, run.runnumber, "tor.dat");
            export_data(host_global.vrad, host_global.particle_number, run.folder_out, run.runnumber, "vrad.dat");
            export_data(host_global.vz,   host_global.particle_number, run.folder_out, run.runnumber, "vz.dat");
            export_data(host_global.vtor, host_global.particle_number, run.folder_out, run.runnumber, "vtor.dat");
            export_table(run.folder_out, run.runnumber, "coords.dat", host_global.particle_number,
                host_global.rad, "R [m]",      host_global.z, "Z [m]",      host_global.tor, "T [m]", 
                host_global.vrad, "v_R [m/s]", host_global.vz, "v_Z [m/s]", host_global.vtor, "v_T [m/s]");
        }*/
        
        printf("\nData folder: %s/%s\n\n", run.folder_out, run.runnumber);
        
        //! Free CUDA
        /*for (int i=0; i<6; ++i){
            cudaFree(host_global.coords[i]);
        }
        cudaFree(host_global.coords);
        cudaFree(dev_common.spline_grid[0]);
        cudaFree(dev_common.spline_grid[1]);
        cudaFree(dev_common.spline_grid);

        for (int i=0; i<3; ++i){
            cudaFree(host_common.espline[i]);
            cudaFree(host_common.bspline[i]);
        }
        cudaFree(host_common.espline);
        cudaFree(host_common.bspline);
        

        //! Free RAM    
        free(dev_common.spline_grid[0]);
        free(dev_common.spline_grid[1]);
        free(dev_common.spline_grid);

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
        
        
        //! FREE SERVICE_VAR variables (RAM, cuda)
        free(SERVICE_VAR);  cudaFree(service_var);
        
        printf("Ready.\n\n");
    }
}

void fill_header_file(ShotProp shot, BeamProp beam, RunProp run, double DETECTOR[5]){
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
    export_header("Detector position (R)", "m", DETECTOR[0], run.folder_out, run.runnumber);
    export_header("Detector position (Z)", "m", DETECTOR[1], run.folder_out, run.runnumber);
    export_header("Detector position (T)", "m", DETECTOR[2], run.folder_out, run.runnumber);
    export_header("Detector angle (Z/R)", "°", atan(DETECTOR[3])*180.0/PI, run.folder_out, run.runnumber);
    export_header("Detector angle (T/R)", "°", atan(DETECTOR[4])*180.0/PI, run.folder_out, run.runnumber);
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
