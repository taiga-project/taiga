#include "feedback.h"
#include "utils/prop.h"
#include "utils/basic_functions.h"
#include "dataio/data_export.h"

void print_run_details(TaigaGlobals *host_global, TaigaCommons *host_common, ShotProp shot, RunProp run){
    printf("===============================\n");
    printf("%s\n", concat("TAIGA ", TAIGA_VERSION," (r", GIT_REV, ")", NULL));
    printf("Shotname: %s\n", shot.name); 
    printf("Detector: %s\n", shot.detector_mask);
    printf("  R:\t%lf\n", host_common->detector_geometry[0]);
    printf("  Z:\t%lf\n", host_common->detector_geometry[1]);
    printf("  T:\t%lf\n", host_common->detector_geometry[2]);
    printf("  angle (Z/R):\t%lf°\n", (host_common->detector_geometry[3])/PI*180.0);
    printf("  angle (T/R):\t%lf°\n", (host_common->detector_geometry[4])/PI*180.0);
    printf("===============================\n");
    printf("Number of blocks (threads): %d\n", run.block_number);
    printf("Block size: %d\n", run.block_size);
    printf("Number of particles: %ld (%ld)\n", run.particle_number, host_global->particle_number);
    printf("Max steps on device (GPU): %ld\n", run.step_device);
    printf("Max steps on host (HDD): %ld\n", run.step_host);
    printf("===============================\n");
}

void fill_header_file(TaigaCommons *common, BeamProp beam, ShotProp shot, RunProp run){
    export_header(concat("TAIGA ", TAIGA_VERSION," (r", GIT_REV, ")", NULL), run.folder_out, run.runnumber);
    export_header_addline(run.folder_out, run.runnumber);
    export_header(concat("Shot ID: ", shot.name, NULL), run.folder_out, run.runnumber);
    export_header(concat("Run ID:  ", run.runnumber, NULL), run.folder_out, run.runnumber);
    export_header_addline(run.folder_out, run.runnumber);
    export_header("ABP ION TRAJECTORIES", run.folder_out, run.runnumber);
    
    if(READINPUTPROF==1){
        export_header("Manual (6D) input profile", run.folder_out, run.runnumber);
    }else if(RENATE==110){
        export_header("TS + Renate 1.1.0 input profile", run.folder_out, run.runnumber);
    }
    export_header_addline(run.folder_out, run.runnumber);
    
    if(!READINPUTPROF){
        export_header(concat("Beam species:\t", beam.species, NULL), run.folder_out, run.runnumber);
        export_header("Beam charge number", "", beam.charge, run.folder_out, run.runnumber);
        export_header("Beam energy", "keV", beam.energy, run.folder_out, run.runnumber);
        export_header("Beam diameter", "mm", beam.diameter*1000, run.folder_out, run.runnumber);
        export_header("Beam deflection (toroidal/vertical)", "°", beam.toroidal_deflection*180.0/PI, beam.vertical_deflection*180.0/PI, run.folder_out, run.runnumber);
    }
    
    export_header("Number of ions", "", run.particle_number, run.folder_out, run.runnumber);
    export_header_addline(run.folder_out, run.runnumber);
    export_header("Detector position (R)", "m", common->detector_geometry[0], run.folder_out, run.runnumber);
    export_header("Detector position (Z)", "m", common->detector_geometry[1], run.folder_out, run.runnumber);
    export_header("Detector position (T)", "m", common->detector_geometry[2], run.folder_out, run.runnumber);
    export_header("Detector angle (Z/R)", "°", (common->detector_geometry[3])*180.0/PI, run.folder_out, run.runnumber);
    export_header("Detector angle (T/R)", "°", (common->detector_geometry[4])*180.0/PI, run.folder_out, run.runnumber);
    export_header(concat("Detector mask:  \t", shot.detector_mask, NULL), run.folder_out, run.runnumber);
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