// TAIGA default parameters

#define $R_defl 2.3                 //! radial position of deflection plates in meter -> TOROIDAL DEFLECTION

#define PI 3.141592653589792346
#define ELEMENTARY_CHARGE 1.60217656535e-19
#define AMU 1.66053892173e-27
#define INFINITY RAND_MAX
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
#include "test/cuda/nvToolsExt.h"

#include "prop.h"
#include "main.cuh"
#include "debug_functions.c"
#include "basic_functions.c"

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

void input_init_taiga(int argc, char *argv[], shot_prop *shot, beam_prop *beam, run_prop *run){
    
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
        }else if (!strcmp(input, "--runnumber_file") || !strcmp(input, "--runnumber") || !strcmp(input, "-r")){
            input = strtok(NULL, "=");
            int runnumber = atoi(input);
            if (runnumber || !strcmp(input, "0")){
                run->runnumber = runnumber;
                strcpy(run->runnumber_file, "console init");
                printf("Runnumber: %d\n", run->runnumber);            
            }else{
                strcpy(run->runnumber_file, input);
                printf("Runnumber file: %s\n", run->runnumber_file);
            }
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
        printf("  -r  --runnumber_file=PATH   Runnumber file path\n");
        printf("  -r  --runnumber=INTEGER     Runnumber value\n");
        printf("  -s, --ion-source=PATH       Ion source path\n");
        printf("      --ion-source-coords=XXX Order of coordinates (RZT or RTZ) in input file\n");
}

int main(int argc, char *argv[]){    
    shot_prop shot; init_shot_prop(&shot);
    beam_prop beam; init_beam_prop(&beam);
    run_prop run;   init_run_prop(&run);
    input_init_taiga(argc, argv, &shot, &beam, &run);
    
    if (run.help == 1){
        print_help_message();
    }else if (run.help == 2){
        set_cuda(1);
    }else{
        device_global dev_global, host_global;
        device_shared dev_shared, host_shared;
        init_device(&dev_global, &dev_shared);
        
        parameter_reader(&shot, &beam, &run);
        runnumber_reader(&shot, &run);
        
        char* folder_out=concat("results/", shot.name);
        char timestamp[80];
        sprintf(timestamp, "%d", run.runnumber);
        
        init_dir(folder_out, timestamp);
        CopyFile(run.parameter_file, concat(folder_out,"/",timestamp,"/parameters.sh"));

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
        
        host_global.particle_number = run.block_size * run.block_number;
        
        if (READINPUTPROF == 1){
            double *X_temp;
            host_global.particle_number = read_vector(&X_temp, "input", "manual_profile", "rad.dat");
            run.block_number = host_global.particle_number / run.block_size+1;
            free(X_temp);
        }
        
        //! CUDA profiler START
        cudaProfilerStart();
        
        set_cuda(run.debug);
        
        // set timestamp
        time_t rawtime;
        struct tm *info;
        
        // coords
        double *X_PTR[3], **x_ptr;
        double *V_PTR[3], **v_ptr;
        size_t dimXP = 3*sizeof(double*);
        
        double *XR,  *xr; 
        double *XZ,  *xz;
        double *XT,  *xt;
        
        double *VR,  *vr; 
        double *VZ,  *vz;
        double *VT,  *vt;
        
        printf("Number of blocks (threads): %d\n", run.block_number);
        printf("Block size: %d\n", run.block_size);
        printf("Number of particles: %d\n", host_global.particle_number);
        printf("Max steps on device (GPU): %d\n", run.step_device);
        printf("Max steps on host (HDD): %d\n", run.step_host);
        
        // phys. constants
        double eperm = ELEMENTARY_CHARGE/ AMU/ beam.mass;
        
       //! coordinates
        init_coords(&host_global, &dev_global, &dev_shared, beam, shot, run);
        
        //! grid
        init_grid(shot, run, &host_shared, &dev_shared);
        
        int magnetic_field_loaded = magnetic_field_read_and_init(shot, run, &host_shared, &dev_shared);
        if (shot.electric_field_module) shot.electric_field_module = electric_field_read_and_init(shot, run, &host_shared, &dev_shared);
        
        // detector cell id
        size_t dimRint = host_global.particle_number * sizeof(int);
        int *DETCELLID, *detcellid;
        DETCELLID = (int *)malloc(dimRint); cudaMalloc((void **) &detcellid,  dimRint);
        
        // service value
        size_t dimService = SERVICE_VAR_LENGTH * sizeof(double);
        double *SERVICE_VAR, *service_var;
        SERVICE_VAR = (double *)malloc(dimService); 
        
        for(int i=0 ; i<SERVICE_VAR_LENGTH ; ++i){
            SERVICE_VAR[i] = 0;
        }
        
        cudaMalloc((void **) &service_var,  dimService);
        cudaMemcpy(service_var, SERVICE_VAR, dimService, cudaMemcpyHostToDevice);
        
printf("212 \n");        //! MEMCOPY (HOST2device)
        
        //! DETECTOR COORDS (HOST2device)
        cudaMemcpy(&(dev_shared.detector_geometry), DETECTOR, dimD, cudaMemcpyHostToDevice);
        
        if (!FASTMODE){
            // OUTPUT INIT
            export_data(host_global.rad,  host_global.particle_number, folder_out, timestamp, "t_rad.dat");
            export_data(host_global.z,    host_global.particle_number, folder_out, timestamp, "t_z.dat");
            export_data(host_global.tor,  host_global.particle_number, folder_out, timestamp, "t_tor.dat");
            export_data(host_global.vrad, host_global.particle_number, folder_out, timestamp, "t_vrad.dat");
            export_data(host_global.vz,   host_global.particle_number, folder_out, timestamp, "t_vz.dat");
            export_data(host_global.vtor, host_global.particle_number, folder_out, timestamp, "t_vtor.dat");
        }
printf("226 ");        
        //! Set CUDA timer 
        cudaEvent_t cuda_event_core_start, cuda_event_core_end, cuda_event_copy_start, cuda_event_copy_end;
        clock_t cpu_event_copy_start, cpu_event_copy_end;
        double cpu_time_copy=0, cuda_time_core, cuda_time_copy;
        float cuda_event_core, cuda_event_copy;
        cudaEventCreate(&cuda_event_core_start);
        cudaEventCreate(&cuda_event_core_end);
        cudaEventCreate(&cuda_event_copy_start);
        cudaEventCreate(&cuda_event_copy_end);
        
        if (run.debug == 1 && !FASTMODE)   debug_message_init(host_global.rad, host_global.z, host_global.tor, host_global.vrad, host_global.vz, host_global.vtor);
        size_t dimX = host_global.particle_number*sizeof(double);
        for (int step_i=0;step_i<run.step_host;step_i++){
printf("257 ");            
            if (step_i == 0) cudaEventRecord(cuda_event_core_start, 0);
           
            taiga <<< run.block_number, run.block_size >>> (dev_global, dev_shared, service_var);
            
            if (step_i == 0) cudaEventRecord(cuda_event_core_end, 0);
            cudaEventSynchronize(cuda_event_core_end);
            //ERRORCHECK();
            
            if (!FASTMODE){
                // ION COORDS (device2HOST)
                if (step_i == 0) cudaEventRecord(cuda_event_copy_start, 0);
                coord_memcopy(&host_global, &dev_global, &dev_shared, beam, shot, run);
                //cudaMemcpy(host_global.rad, device_global.rad, dimX, cudaMemcpyDeviceToHost);
               /* cudaMemcpy(XZ, xz, dimX, cudaMemcpyDeviceToHost);
                cudaMemcpy(XT, xt, dimX, cudaMemcpyDeviceToHost);
                cudaMemcpy(VR, vr, dimX, cudaMemcpyDeviceToHost);
                cudaMemcpy(VZ, vz, dimX, cudaMemcpyDeviceToHost);
                cudaMemcpy(VT, vt, dimX, cudaMemcpyDeviceToHost);*/
                //ERRORCHECK();                
                if (step_i == 0) cudaEventRecord(cuda_event_copy_end, 0);
                
printf("278 ");                // Save data to files
                cpu_event_copy_start = clock();
                export_data(host_global.rad,  host_global.particle_number, folder_out, timestamp, "t_rad.dat");
                export_data(host_global.z,    host_global.particle_number, folder_out, timestamp, "t_z.dat");
                export_data(host_global.tor,  host_global.particle_number, folder_out, timestamp, "t_tor.dat");
                export_data(host_global.vrad, host_global.particle_number, folder_out, timestamp, "t_vrad.dat");
                export_data(host_global.vz,   host_global.particle_number, folder_out, timestamp, "t_vz.dat");
                export_data(host_global.vtor, host_global.particle_number, folder_out, timestamp, "t_vtor.dat");
                cpu_event_copy_end = clock();
            }
printf("288 ");            
            if (run.debug == 1)    printf("Step\t%d/%d\n",step_i,run.step_host);
            if (run.debug == 1 && !FASTMODE)    debug_message_run(host_global.rad, host_global.z, host_global.tor, host_global.vrad, host_global.vz, host_global.vtor);
        }
        
        // Get CUDA timer 
        cudaEventElapsedTime(&cuda_event_core, cuda_event_core_start, cuda_event_core_end);
        cudaEventElapsedTime(&cuda_event_copy, cuda_event_copy_start, cuda_event_copy_end);
        if (!FASTMODE) cpu_time_copy = ((double) (4.0+run.step_host)*(cpu_event_copy_end - cpu_event_copy_start)) / CLOCKS_PER_SEC;
        cuda_time_copy = (double) (1.0+run.step_host)*cuda_event_copy/1000.0;
        cuda_time_core =  run.step_host*cuda_event_core/1000.0;
        
        printf("===============================\n");
        printf ("CUDA kernel runtime: %lf s\n", cuda_time_core);
        printf ("CUDA memcopy time:   %lf s\n", cuda_time_copy);
        printf ("CPU->HDD copy time:  %lf s\n", cpu_time_copy);    
        printf("===============================\n");
        
        undetected <<<1,1>>>(detcellid, host_global.particle_number, service_var);
        
        //! MEMCOPY (device2HOST)
        cudaMemcpy(SERVICE_VAR, service_var, dimService, cudaMemcpyDeviceToHost);
        if(SERVICE_VAR[0]!=42.24){
            printf("\n +----------------------------+\n | Fatal error in running.    | \n | The CUDA did not run well. |\n | Service value: %11lf |\n +----------------------------+\n\n", SERVICE_VAR[0]);
        }else{
            printf("\nSuccessful run. \n\n");
        }
        
        printf("Lost particle ratio: \t %.4lf % \n", SERVICE_VAR[1]*100);
        
        detector_module(x_ptr, detector, detcellid, shot.detector_mask, run.block_number, run.block_size, host_global.particle_number, folder_out, timestamp);
        cudaMemcpy(DETCELLID, detcellid, dimRint, cudaMemcpyDeviceToHost);
        export_data(DETCELLID, host_global.particle_number, folder_out, timestamp, "detector", "cellid.dat");
        
        if (run.debug == 1)    debug_service_vars(SERVICE_VAR);
        
        //! CUDA profiler STOP
        cudaProfilerStop();
        
        export_header(concat("TAIGA ", TAIGA_VERSION," (r", GIT_REV, ")"), folder_out, timestamp);
        export_header_addline(folder_out, timestamp);
        export_header(concat("Shot ID: ",shot.name), folder_out, timestamp);
        export_header(concat("Run ID:  ",timestamp), folder_out, timestamp);
        export_header_addline(folder_out, timestamp);
        export_header("ABP ION TRAJECTORIES", folder_out, timestamp);
        
        if(READINPUTPROF==1){
            export_header("Manual (6D) input profile", folder_out, timestamp);
        }else if(RENATE==110){
            export_header("TS + Renate 1.1.0 input profile", folder_out, timestamp);
        }
        export_header_addline(folder_out, timestamp);
        
        if(!READINPUTPROF){
            export_header("Beam energy", "keV", beam.energy, folder_out, timestamp);
            export_header("Atomic mass", "AMU", beam.mass, folder_out, timestamp);
            export_header("Beam diameter", "mm", beam.diameter*1000, folder_out, timestamp);
            export_header("Beam deflection (toroidal/vertical)", "°", beam.toroidal_deflection*180.0/PI, beam.vertical_deflection*180.0/PI, folder_out, timestamp);
        }
        
        export_header("Number of ions", "", (double)host_global.particle_number, folder_out, timestamp);
        export_header_addline(folder_out, timestamp);
        export_header("Detector position (R)", "m", DETECTOR[0], folder_out, timestamp);
        export_header("Detector position (Z)", "m", DETECTOR[1], folder_out, timestamp);
        export_header("Detector position (T)", "m", DETECTOR[2], folder_out, timestamp);
        export_header("Detector angle (Z/R)", "°", atan(DETECTOR[3])*180.0/PI, folder_out, timestamp);
        export_header("Detector angle (T/R)", "°", atan(DETECTOR[4])*180.0/PI, folder_out, timestamp);
        export_header(concat("Detector mask:  \t", shot.detector_mask), folder_out, timestamp);
        export_header_addline(folder_out, timestamp);
        export_header("Timestep", "s", run.timestep, folder_out, timestamp);
        export_header_addline(folder_out, timestamp);
        export_header("CUDA kernel runtime", "s", cuda_time_core, folder_out, timestamp);
        export_header("CUDA memcopy time", "s", cuda_time_copy, folder_out, timestamp);
        export_header("CPU->HDD copy time", "s", cpu_time_copy, folder_out, timestamp);
        export_header_addline(folder_out, timestamp);
        export_header("Number of blocks (threads)", "", run.block_number, folder_out, timestamp);
        export_header("Block size", "", run.block_size, folder_out, timestamp);
        export_header("Length of a loop", "", run.step_device, folder_out, timestamp);
        export_header("Number of loops", "", run.step_host, folder_out, timestamp);
        
        if (!FASTMODE){
            //! Save data to files
            export_data(host_global.rad, host_global.particle_number, folder_out, timestamp, "rad.dat");
            export_data(host_global.z,   host_global.particle_number, folder_out, timestamp, "z.dat");
            export_data(host_global.tor, host_global.particle_number, folder_out, timestamp, "tor.dat");
            export_data(host_global.vrad, host_global.particle_number, folder_out, timestamp, "vrad.dat");
            export_data(host_global.vz,   host_global.particle_number, folder_out, timestamp, "vz.dat");
            export_data(host_global.vtor, host_global.particle_number, folder_out, timestamp, "vtor.dat");
            export_table(folder_out, timestamp, "coords.dat", host_global.particle_number,
                host_global.rad, "R [m]",      host_global.z, "Z [m]",      host_global.tor, "T [m]", 
                host_global.vrad, "v_R [m/s]", host_global.vz, "v_Z [m/s]", host_global.vtor, "v_T [m/s]");
        }
        
        printf("\n\nData folder: %s/%s\n\n", folder_out, timestamp);
        
        //! Free CUDA
        /*#for (int i=0; i<6; ++i){
            cudaFree(host_global.coords[i]);
        }
        cudaFree(host_global.coords);
        cudaFree(dev_shared.spline_grid[0]);
        cudaFree(dev_shared.spline_grid[1]);
        cudaFree(dev_shared.spline_grid);

        for (int i=0; i<3; ++i){
            cudaFree(host_shared.espline[i]);
            cudaFree(host_shared.bspline[i]);
        }
        cudaFree(host_shared.espline);
        cudaFree(host_shared.bspline);
        

        //! Free RAM    
        free(dev_shared.spline_grid[0]);
        free(dev_shared.spline_grid[1]);
        free(dev_shared.spline_grid);

        for (int i=0; i<3; ++i){
            free(host_shared.espline[i]);
            free(host_shared.bspline[i]);
        }
        free(host_shared.espline);
        free(host_shared.bspline);
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
