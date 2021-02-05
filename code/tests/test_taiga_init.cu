#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
//#include <stdbool.h>
//#include <curand_kernel.h>
#include "../basic_functions.h"
#include "../taiga_constants.h"
#include "../prop.h"
#include "../taiga_constants.h"
#include "../running/generate_coords.cu"
#include "../dataio/data_import.c"
#include "../taiga_init.c"
#include "../dataio/parameter_reader.c"

#include "../dataio/beam.h"
#if READINPUTPROF == 1
    #include "../dataio/beam_manual_profile.c"
#elif RENATE == 110
    #include "../dataio/beam_renate110.c"
#else
    #error A valid beam module is required!
#endif

#include "../running/detector_postproc.cu"
#include "../running/detector_sum.cu"

#define LENGTH_TMP 10

__global__ void test_init_grid_cuda(TaigaCommons* c, double *tmp){
    tmp[0] = PI;
    tmp[1] = c->grid_size[0];
    tmp[2] = c->grid_size[1];
    tmp[3] = c->spline_rgrid[0];
    tmp[4] = c->spline_rgrid[c->grid_size[0]-1];
    tmp[5] = c->spline_zgrid[0];
    tmp[6] = c->spline_zgrid[c->grid_size[1]-1];
    /*tmp[7] = c->detector_geometry[0];
    tmp[8] = c->detector_geometry[1];
    tmp[9] = c->detector_geometry[2];*/
}

__global__ void test_init_coords_cuda(TaigaGlobals* g, double *tmp){
    tmp[0] = g->particle_number;
    tmp[1] = g->rad[0];
    tmp[2] = g->z[0];
    tmp[3] = g->tor[0];
    tmp[4] = g->rad[1];
    tmp[5] = g->z[1];
    tmp[6] = g->tor[1];
    tmp[7] = g->rad[g->particle_number-1];
    tmp[8] = g->z[g->particle_number-1];
    tmp[9] = g->tor[g->particle_number-1];
}

void init_tmp(double *h_tmp){
    for (int i=0; i<LENGTH_TMP; ++i){
        h_tmp[i] = UNDEFINED_FLOAT;
    }
}

void print_tmp(double *h_tmp){
    for (int i=0; i<LENGTH_TMP; ++i){
        printf("TMP %d: %lf\n", i, h_tmp[i]);
    }
}

void start_reference(double **h_tmp, double **d_tmp){
    size_t dim_tmp = sizeof(double)*LENGTH_TMP;
    *h_tmp = (double *) malloc(dim_tmp);
    init_tmp(*h_tmp);
    cudaMalloc((void **) d_tmp, dim_tmp);
    cudaMemcpy(*d_tmp, *h_tmp, dim_tmp, cudaMemcpyHostToDevice);
}

void end_reference(double **h_tmp, double **d_tmp){
    size_t dim_tmp = sizeof(double)*LENGTH_TMP;
    cudaMemcpy(*h_tmp, *d_tmp, dim_tmp, cudaMemcpyDeviceToHost);
    print_tmp(*h_tmp);
}

void test_init_grid(){
    printf("Test: test_init_grid()l()\n");
    int tmp_length = 20;

    ShotProp shot;
    RunProp run;
    TaigaCommons *host_common, *shared_common, *dev_common;
    
    size_t dim_commons = sizeof(TaigaCommons);
    host_common = (TaigaCommons*)malloc(dim_commons);
    shared_common = (TaigaCommons*)malloc(dim_commons);
    cudaMalloc((void **) &dev_common, dim_commons);
    
    strcpy(shot.name, "17178_1097");
    printf("Init grid\n");
    init_grid(shot, run, host_common, shared_common);
    cudaMemcpy(dev_common, shared_common, dim_commons, cudaMemcpyHostToDevice);
        
    double *h_tmp, *d_tmp;
    start_reference(&h_tmp, &d_tmp);
    test_init_grid_cuda <<< 1, 1 >>> (dev_common, d_tmp);
    end_reference(&h_tmp, &d_tmp);
    
    printf("R = %lf ... %lf\n", host_common->spline_rgrid[0], host_common->spline_rgrid[host_common->grid_size[0]-1]);
    printf("Z = %lf ... %lf\n", host_common->spline_zgrid[0], host_common->spline_zgrid[host_common->grid_size[1]-1]);
}

void test_init_coords(){
    printf("Test: test_init_coords()\n");
    
    ShotProp shot; init_shot_prop(&shot);
    BeamProp beam; init_beam_prop(&beam);
    RunProp run;   init_run_prop(&run);
    TaigaGlobals *host_global, *shared_global, *dev_global;
    
    size_t size_globals = sizeof(TaigaGlobals);
    host_global = (TaigaGlobals*)malloc(size_globals);
    shared_global = (TaigaGlobals*)malloc(size_globals);
    cudaMalloc((void **) &dev_global, size_globals);
    
    strcpy(shot.name, "17178_1097");
    init_taiga_props("particles", "1000", &beam, &shot, &run);
    
    init_coords(&beam, &shot, &run, host_global, shared_global);
    set_particle_number(&run, host_global, shared_global);
    cudaMemcpy(dev_global, shared_global, size_globals, cudaMemcpyHostToDevice);
        
    double *h_tmp, *d_tmp;
    start_reference(&h_tmp, &d_tmp);
    test_init_coords_cuda <<< 1, 1 >>> (dev_global, d_tmp);
    end_reference(&h_tmp, &d_tmp);
    
    printf("Number of ions: %ld\n", host_global->particle_number);
    
    printf("1st:  [%lf, %lf, %lf]\n", host_global->rad[0], host_global->z[0], host_global->tor[0]);
    printf("2nd:  [%lf, %lf, %lf]\n", host_global->rad[1], host_global->z[1], host_global->tor[1]);
    printf("Last: [%lf, %lf, %lf]\n", host_global->rad[host_global->particle_number-1], host_global->z[host_global->particle_number-1], host_global->tor[host_global->particle_number-1]);
}

__global__ void test_detector_cuda(DetectorProp *d,  double *tmp){
    tmp[0] = PI;
    tmp[1] = d->xgrid[0];
    tmp[2] = d->ygrid[0];
}

void test_init_detector(){
    printf("Test: test_detector()\n");
    ShotProp shot; init_shot_prop(&shot);
    DetectorProp *shared_detector, *device_detector;
    
    size_t size_detector_prop = sizeof(DetectorProp);
    shared_detector = (DetectorProp*)malloc(size_detector_prop);
    cudaMalloc((void **) &device_detector, size_detector_prop);
    
    init_detector(shared_detector, device_detector, shot);
    
    double *h_tmp, *d_tmp;
    start_reference(&h_tmp, &d_tmp);
    test_detector_cuda <<< 1, 1 >>> (device_detector, d_tmp);
    end_reference(&h_tmp, &d_tmp);
}

__global__ void test_detector_struct(TaigaCommons *c, DetectorProp *d, double *tmp){
    tmp[0] = c->detector_geometry[0];
    tmp[1] = c->detector_geometry[1];
    tmp[2] = c->detector_geometry[2];
    tmp[3] = c->detector_geometry[3];
    tmp[4] = c->detector_geometry[4];
    tmp[5] = d->xgrid[0];
    tmp[6] = d->ygrid[0];
    tmp[7] = d->number_of_detector_cells;
}

__global__ void test_detector_detcellid(TaigaGlobals *g, double *tmp){
    tmp[0] = g->detcellid[0];
    tmp[1] = g->detcellid[1];
    tmp[2] = g->detcellid[2];
    tmp[3] = g->detcellid[3];
    tmp[4] = g->detcellid[4];
    tmp[5] = g->detcellid[5];
    tmp[6] = g->particle_number;
}

__global__ void test_detector_full_cuda(DetectorProp *d, double *tmp){
    tmp[0] = d->counter[0];
    tmp[1] = d->counter[1];
    tmp[2] = d->counter[2];
    tmp[3] = d->counter[9];
    tmp[4] = d->counter[10];
    tmp[5] = d->counter[11];
    tmp[6] = d->counter[12];
    tmp[7] = d->counter[13];
    tmp[8] = d->counter[14];
    tmp[9] = d->counter[15];
}

void test_init_detector_full(){
    printf("Test: test_init_detector_full()\n");
    double z_unit = 0.001;
    double z_center = 0.2101;
    double tor_unit = 0.001;
    double tor_center = 0.0101;

    ShotProp shot; init_shot_prop(&shot);
    RunProp run;   init_run_prop(&run);
    TaigaCommons *host_common, *shared_common, *device_common;
    TaigaGlobals *host_global, *shared_global, *device_global;
    DetectorProp *shared_detector, *device_detector;
    
    size_t dim_commons = sizeof(TaigaCommons);
    host_common = (TaigaCommons*)malloc(dim_commons);
    shared_common = (TaigaCommons*)malloc(dim_commons);
    cudaMalloc((void **) &device_common, dim_commons);
    
    strcpy(shot.name, "17178_1097");
    strcpy(shot.detector_geometry, "0.7, 0.2, 0, 0, 0");
    run.block_number = 20;
    run.block_size = 20;
    
    init_grid(shot, run, host_common, shared_common);
    set_detector_geometry(shot, host_common, shared_common);
    printf("Detector geometry: [%lf %lf %lf %lf %lf]\n",
        host_common->detector_geometry[0],
        host_common->detector_geometry[1],
        host_common->detector_geometry[2],
        host_common->detector_geometry[3],
        host_common->detector_geometry[4]);
    
    
    size_t size_detector_prop = sizeof(DetectorProp);
    shared_detector = (DetectorProp*)malloc(size_detector_prop);
    cudaMalloc((void **) &device_detector, size_detector_prop);
    
    init_detector(shared_detector, device_detector, shot);
    
    size_t size_globals = sizeof(TaigaGlobals);
    host_global = (TaigaGlobals*)malloc(size_globals);
    shared_global = (TaigaGlobals*)malloc(size_globals);
    cudaMalloc((void **) &device_global, size_globals);
    
    size_t size_coord = run.block_size * run.block_number * sizeof(double);
    size_t size_detcellid = run.block_size * run.block_number * sizeof(int);
    host_global->rad = (double*)malloc(size_coord);
    host_global->z   = (double*)malloc(size_coord);
    host_global->tor = (double*)malloc(size_coord);
    host_global->detcellid = (int*)malloc(size_detcellid);
//    host_global->vrad = (double*)malloc(size_coord);
//    host_global->vz   = (double*)malloc(size_coord);
//    host_global->vtor = (double*)malloc(size_coord);
    
    for (int i=0; i<run.block_number; ++i){
        for (int j=0; j<run.block_size; ++j){
            int k = i*run.block_size + j;
            host_global->rad[k] = 0.7;
            host_global->z[k]   = (i-run.block_number/2)*z_unit+z_center;
            host_global->tor[k] = (j-run.block_size/2)*tor_unit+tor_center;
            //printf("particles: [%lf %lf]\n", host_global->tor[k]*1000, host_global->z[k]*1000);
            host_global->detcellid[k] = 0;
        }
    }
    host_global->particle_number = run.block_size*run.block_number;
    double* shared_rad;
    double* shared_z;
    double* shared_tor;
//    double* shared_vrad;
//    double* shared_vz;
//    double* shared_vtor;
    int* shared_detcellid;
    
    memcpy(shared_global, host_global, size_globals);
    
    cudaMalloc((void **) &shared_rad, size_coord);
    cudaMalloc((void **) &shared_z,   size_coord);
    cudaMalloc((void **) &shared_tor, size_coord);
//    cudaMalloc((void **) &shared_vrad, size_coord);
//    cudaMalloc((void **) &shared_vz,   size_coord);
//    cudaMalloc((void **) &shared_vtor, size_coord);
    cudaMalloc((void **) &shared_detcellid, size_detcellid);
    
    cudaMemcpy(shared_rad,       host_global->rad,       size_coord,  cudaMemcpyHostToDevice);
    cudaMemcpy(shared_z,         host_global->z,         size_coord,  cudaMemcpyHostToDevice);
    cudaMemcpy(shared_tor,       host_global->tor,       size_coord,  cudaMemcpyHostToDevice);
//    cudaMemcpy(shared_vrad,      host_global->vrad,      size_coord,  cudaMemcpyHostToDevice);
//    cudaMemcpy(shared_vz,        host_global->vz,        size_coord,  cudaMemcpyHostToDevice);
//    cudaMemcpy(shared_vtor,      host_global->vtor,      size_coord,  cudaMemcpyHostToDevice);
    cudaMemcpy(shared_detcellid, host_global->detcellid, size_detcellid, cudaMemcpyHostToDevice);
    
    shared_global->rad  = shared_rad;
    shared_global->z    = shared_z;
    shared_global->tor  = shared_tor;
//    shared_global->vrad = shared_vrad;
//    shared_global->vz   = shared_vz;
//    shared_global->vtor = shared_vtor;
    shared_global->detcellid = shared_detcellid;
    cudaMemcpy(device_global, shared_global, sizeof(TaigaGlobals), cudaMemcpyHostToDevice);
    cudaMemcpy(device_common, shared_common, sizeof(TaigaCommons), cudaMemcpyHostToDevice);
    
    detector_postproc <<< run.block_number, run.block_size >>> (device_global, device_common, device_detector);
//    detector_sum <<<1,1>>> (device_global, device_common, device_detector);
    
    double *h_tmp, *d_tmp;
    start_reference(&h_tmp, &d_tmp);
    
    printf("test_detector_detcellid\n");
    test_detector_detcellid <<< 1, 1 >>> (device_global, d_tmp);
    end_reference(&h_tmp, &d_tmp);
    
    printf("test_detector_struct\n");
    test_detector_struct <<< 1, 1 >>> (device_common, device_detector, d_tmp);
    end_reference(&h_tmp, &d_tmp);
    
    printf("test_detector_full_cuda\n");
    test_detector_full_cuda <<< 1, 1 >>> (device_detector, d_tmp);
    end_reference(&h_tmp, &d_tmp);
}

__global__ void test_detector_range_cuda(double *tmp){
    
    tmp[7] = is_in_range(tmp[4], tmp[0], tmp[3]);
    tmp[8] = is_in_range(tmp[5], tmp[0], tmp[3]);
    tmp[9] = is_in_range(tmp[6], tmp[0], tmp[3]);
}

__global__ void test_detector_index_cuda(double *tmp){
    
    double array[] = {tmp[0], tmp[1], tmp[1], tmp[2], tmp[2], tmp[3]};
    
    tmp[7] = get_cell_array_index(tmp[4], array, 3);
    tmp[8] = get_cell_array_index(tmp[5], array, 3);
    tmp[9] = get_cell_array_index(tmp[6], array, 3);
}

void test_detector_conversion(){
    printf("Test: test_detector_conversion()\n");
    double *h_tmp, *d_tmp;
    size_t dim_tmp = sizeof(double)*LENGTH_TMP;
    h_tmp = (double *) malloc(dim_tmp);
    init_tmp(h_tmp);
    
    h_tmp[0] = 0.20;
    h_tmp[1] = 0.21;
    h_tmp[2] = 0.22;
    h_tmp[3] = 0.23;
    h_tmp[4] = 0.195;
    h_tmp[5] = 0.224;
    h_tmp[6] = 0.242;
    
    cudaMalloc((void **) &(d_tmp), dim_tmp);
    cudaMemcpy(d_tmp, h_tmp, dim_tmp, cudaMemcpyHostToDevice);
    
    test_detector_range_cuda <<< 1, 1 >>> (d_tmp);
    end_reference(&h_tmp, &d_tmp);
    
    test_detector_index_cuda <<< 1, 1 >>> (d_tmp);
    end_reference(&h_tmp, &d_tmp);
}
