#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "utils/taiga_constants.h"
#include "utils/prop.c"
#include "main.cuh"
#include "interface/data_export/save_data.c"
#include "interface/feedback.c"
#include "utils/debug_functions.c"
#include "utils/basic_functions.h"
#include "utils/dataio/dir_functions.c"
#include "utils/cuda.cu"
#include "utils/physics.c"
#include "utils/free.cu"

#include "utils/dataio/data_import.c"
#include "interface/data_import/field_import.cu"
#include "interface/parameter_reader.c"
#include "interface/compiler_definition.c"

#include "init/structures/beam.cu"
#include "init/device/init.cu"
#include "init/device/sync.cu"
#include "init/structures/detector.cu"
#include "init/fast_mode.cu"
#include "init/thomson.cu"

#include "interface/data_import/beam_manual_profile.c"
#include "interface/data_import/beam_renate.c"

#include "utils/dataio/data_export.c"

#include "core/maths/maths.cu"
#include "core/solvers/rk4.cu"
#include "core/solvers/runge_kutta_nystrom.cu"
#include "core/solvers/solvers.cuh"
#include "core/solvers/verlet.cu"
#include "core/solvers/boris.cu"
#include "core/solvers/yoshida.cu"
#include "core/detection.cu"
#include "core/maths/cyl2tor.cu"
#include "core/localise_field.cu"
#include "core/maths/bspline.cu"
#include "core/traj.cu"
#include "core/init/generate_coords.cu"
#include "core/taiga.cu"
#include "core/init/init_beamlet.cu"
#include "core/physics/ionisation.cu"

#include "detector/module.cu"
#include "detector/postproc.cu"
#include "detector/sum.cu"

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
            run->help = HELP_MODE;
        }else if (!strcmp(input, "--devices") || !strcmp(input, "-D") || !strcmp(input, "-l")){
            run->help = HELP_DEVICES;
        }else if (!strcmp(input, "--parameter_file") || !strcmp(input, "-p")){
            input = strtok(NULL, "=");
            strcpy(run->parameter_file, input);
            printf("Parameter file: %s\n", run->parameter_file);
        }else if (!strcmp(input, "--runnumber_file")  || !strcmp(input, "-R")){
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
        }else if (!strcmp(input, "--ion-source-coords") || !strcmp(input, "-S")){
            input = strtok(NULL, "=");
            strcpy(run->io_coordinate_order, input);
            printf("Order of coordinates in input file: %s\n", run->io_coordinate_order);
        }else if (!strcmp(input, "--version") || !strcmp(input, "-v")){
            input = strtok(NULL, "=");
            run->help = HELP_VERSION;
        }else{
            printf("Warning: Undefined command line parameter: %s\n", input);
        }
    }
}

void print_help_message(){
    printf("%s\n", concat("TAIGA ", TAIGA_VERSION," (r", GIT_REV, ")", NULL));
    printf("Usage: taiga.exe [OPTION]\nOptions:\n");
    printf("  -d,      --debug                 Print additional debug informations\n");
    printf("  -D, -l,  --devices               List GPU devices\n");
    printf("  -f,      --fulltrace             Save coordinates at every timestep\n");
    printf("  -h,      --help                  Help message\n");
    printf("  -p=PATH, --parameter_file=PATH   Parameter file path\n");
    printf("  -r=INT,  --runnumber=INTEGER     Runnumber value\n");
    printf("  -R=PATH  --runnumber_file=PATH   Runnumber file path\n");
    printf("  -s=PATH, --ion-source=PATH       Ion source path\n");
    printf("  -S=XXX   --ion-source-coords=XXX Order of coordinates (RZT or RTZ) in input file\n");
    printf("  -v       --version               Version number\n");
}

void print_version(){
    printf("TAIGA (%s)\n\n", TAIGA_VERSION);
    printf("Trajectory simulator of ABP Ions with GPU Acceleration\n");
    printf("Copyright (C) 2011--2024\n\n");
    printf("Written by Matyas Aradi\n");
}

int main(int argc, char *argv[]){
    ShotProp shot; init_shot_prop(&shot);
    BeamProp beam; init_beam_prop(&beam);
    RunProp run;   init_run_prop(&run);
    input_init_taiga(argc, argv, &shot, &beam, &run);
    
    if (run.help == HELP_MODE){
        print_help_message();
    }else if (run.help == HELP_DEVICES){
        set_cuda(1);
    }else if (run.help == HELP_VERSION){
        print_version();
    }else{
        read_compiler_definition(&run);
        parameter_reader(&beam, &shot, &run);
        runnumber_reader(&shot, &run);
        
        init_dir(run.folder_out, run.runnumber);
        CopyFile(run.parameter_file, concat(run.folder_out,"/",run.runnumber,"/parameters.sh", NULL));
        
        set_cuda(run.debug);

        //! Set CUDA timer
        cudaEvent_t cuda_event_core_start, cuda_event_core_end, cuda_event_copy_start, cuda_event_copy_end;
        clock_t cpu_event_copy_start, cpu_event_copy_end;
        float cuda_event_core, cuda_event_copy;
        cudaEventCreate(&cuda_event_core_start);
        cudaEventCreate(&cuda_event_core_end);
        cudaEventCreate(&cuda_event_copy_start);
        cudaEventCreate(&cuda_event_copy_end);

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

        CHECK_ERROR(cudaMalloc((void **) &device_global, size_global));
        CHECK_ERROR(cudaMalloc((void **) &device_common, size_commons));
        CHECK_ERROR(cudaMalloc((void **) &device_detector, size_detector_prop));
        
        init_host(host_global, host_common);
        
        set_particle_number(&run, host_global, shared_global);
        
        //! coordinates
        init_coords(&beam, &shot, &run, host_global, shared_global);
        
        //! grid
        init_grid(shot, run, host_common, shared_common);
        switch (run.field_interpolation_method) {
            case CUBIC_SPLINE:
                magnetic_field_read_and_init(shot, run, host_common, shared_common);
                poloidal_flux_read_and_init(shot, run, host_common, shared_common);
                break;
            case CUBIC_BSPLINE:
                magnetic_field_read_and_init_with_bsplines(shot, run, host_common, shared_common);
                poloidal_flux_read_and_init_with_bsplines(shot, run, host_common, shared_common);
                break;
            default:
                printf("Invalid interpolation value\n");
                exit(1);
        }

        if (run.is_electric_field_on) run.is_electric_field_on = electric_field_read_and_init(shot, run, host_common, shared_common);
        if (run.is_magnetic_field_perturbation) run.is_magnetic_field_perturbation = poloidal_flux_read_and_init(shot, run, host_common, shared_common);

        // detector
        set_detector_geometry(shot, host_common, shared_common, shared_detector);
        if (shared_detector->detector_module_on) {
            init_detector(shared_detector, device_detector, shot);
        }

        // <service value>
        size_t dimService = SERVICE_VAR_LENGTH * sizeof(double);
        double *host_service_array, *device_service_array;
        host_service_array = (double *)malloc(dimService);
        
        for(int i=0; i<SERVICE_VAR_LENGTH; ++i){
            host_service_array[i] = 0;
        }
        
        host_service_array[4] = 55555.55555;
        CHECK_ERROR(cudaMalloc((void **) &device_service_array,  dimService));
        CHECK_ERROR(cudaMemcpy(device_service_array, host_service_array, dimService, cudaMemcpyHostToDevice));
        // </service value>
        
        if (run.mode == ALL_IO){
           save_trajectories(host_global, run);
        }
        
        print_run_details(host_global, host_common, shot, run);

        if (run.debug == 1 && run.mode == ALL_IO)   debug_message_init(host_global);

        set_thomson_profiles(shot, host_common, shared_common);

        init_device_structs(beam, shot, run, shared_global, shared_common);
        sync_device_structs(device_global, shared_global, device_common, shared_common, run.mode == ALL_IO);
        if (run.mode != ALL_IO)   init_fastmode(beam, shot, run, device_global);
        
        for (long step_i=0; step_i<run.step_host; ++step_i){
            if (step_i == 0) cudaEventRecord(cuda_event_core_start, 0);
            
            taiga <<< run.block_number, run.block_size >>> (device_global, device_common, device_service_array);
            
            if (step_i == 0) cudaEventRecord(cuda_event_core_end, 0);
            CHECK_ERROR(cudaEventSynchronize(cuda_event_core_end));
            
            if (run.mode == ALL_IO){
                // ION COORDS (device2HOST)
                if (step_i == 0) cudaEventRecord(cuda_event_copy_start, 0);
                coord_memcopy_back(beam, shot, run, host_global, shared_global);
                if (step_i == 0) cudaEventRecord(cuda_event_copy_end, 0);
                
                // Save data to files
                cpu_event_copy_start = clock();
                save_trajectories(host_global, run);
                cpu_event_copy_end = clock();
            }
            
            if (run.debug == 1)    printf("Step\t%ld/%ld\n",step_i,run.step_host);
            if (run.debug == 1 && run.mode == ALL_IO)    debug_message_run(host_global);
        }
        
        // Get CUDA timer
        cudaEventElapsedTime(&cuda_event_core, cuda_event_core_start, cuda_event_core_end);
        cudaEventElapsedTime(&cuda_event_copy, cuda_event_copy_start, cuda_event_copy_end);
        if (run.mode == ALL_IO) run.cpu_time_copy = ((double) (4+run.step_host)*(double)(cpu_event_copy_end - cpu_event_copy_start)) / CLOCKS_PER_SEC;
        run.cuda_time_copy = (double) (1+run.step_host)*cuda_event_copy/1000.0;
        run.cuda_time_core =  (double)run.step_host*cuda_event_core/1000.0;
        printf("===============================\n");
        printf ("CUDA kernel runtime: %lf s\n", run.cuda_time_core);
        printf ("CUDA memcopy time:   %lf s\n", run.cuda_time_copy);
        if (run.mode == ALL_IO)  printf ("CPU->HDD copy time:  %lf s\n", run.cpu_time_copy);
        printf("===============================\n");
        
        //! MEMCOPY (device2HOST)
        CHECK_ERROR(cudaMemcpy(host_service_array, device_service_array, dimService, cudaMemcpyDeviceToHost));
        if(host_service_array[0] != 42.24){
            printf("\n +----------------------------+\n | Fatal error in running.    | \n | The CUDA did not run well. |\n | Service value: %11lf |\n +----------------------------+\n\n", host_service_array[0]);
        }else{
            printf("\nSuccessful run. \n\n");
        }
        
        detector_postproc <<< run.block_number, run.block_size >>> (device_global, device_common, device_detector);
        detector_sum <<<1,1>>> (device_global, device_common, device_detector);
        export_detector(shared_detector, device_detector, shared_global, shot, run);

        if (run.debug == 1)    debug_service_vars(host_service_array);
        
        fill_header_file(host_common, beam, shot, run);
        
        if (run.mode == ALL_IO){
            save_endpoints(host_global, run);
        }else{
            printf("Warning: End-points are not saved in FASTMODE\n");
        }
        
        printf("\nData folder: %s/%s\n\n", run.folder_out, run.runnumber);
        
        //! FREE
        free_taiga(host_global, shared_global, device_global,
                   host_common, shared_common, device_common,
                   shared_detector, device_detector,
                   &shot, &beam, &run,
                   host_service_array, device_service_array);

        CHECK_ERROR(cudaThreadExit());
        printf("Ready.\n\n");
    }
}
