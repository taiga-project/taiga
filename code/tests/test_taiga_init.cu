#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
//#include <stdbool.h>
//#include <curand_kernel.h>
#include "../basic_functions.h"
#include "../prop.h"
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


#define LENGTH_TMP 10

__global__ void test_init_grid_cuda(TaigaCommons* c, double *tmp){
    tmp[0] = PI;
    tmp[1] = c->grid_size[0];
    tmp[2] = c->grid_size[1];
    tmp[3] = c->spline_rgrid[0];
    tmp[4] = c->spline_rgrid[c->grid_size[0]-1];
    tmp[5] = c->spline_zgrid[0];
    tmp[6] = c->spline_zgrid[c->grid_size[1]-1];
}

__global__ void test_init_coords_cuda(TaigaGlobals* g, double *tmp){
    tmp[0] = g->particle_number;
    tmp[1] = g->rad[0];
    tmp[2] = g->z[0];
    tmp[3] = g->tor[0];
    tmp[4] = PI;
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

void test_init_grid(){
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
    size_t dim_tmp = sizeof(double)*LENGTH_TMP;
    h_tmp = (double *) malloc(dim_tmp);
    init_tmp(h_tmp);
    cudaMalloc((void **) &(d_tmp), dim_tmp);
    cudaMemcpy(d_tmp, h_tmp, dim_tmp, cudaMemcpyHostToDevice);
    
    test_init_grid_cuda <<< 1, 1 >>> (dev_common, d_tmp);
    
    cudaMemcpy(h_tmp, d_tmp, dim_tmp, cudaMemcpyDeviceToHost);
    
    print_tmp(h_tmp);
    
    printf("R = %lf ... %lf\n", host_common->spline_rgrid[0], host_common->spline_rgrid[host_common->grid_size[0]-1]);
    printf("Z = %lf ... %lf\n", host_common->spline_zgrid[0], host_common->spline_zgrid[host_common->grid_size[1]-1]);
}

void test_init_coords(){
    int tmp_length = 20;

    ShotProp shot; init_shot_prop(&shot);
    BeamProp beam; init_beam_prop(&beam);
    RunProp run;   init_run_prop(&run);
    TaigaGlobals *host_global, *shared_global, *dev_global;
    
    size_t dim_globals = sizeof(TaigaGlobals);
    host_global = (TaigaGlobals*)malloc(dim_globals);
    shared_global = (TaigaGlobals*)malloc(dim_globals);
    cudaMalloc((void **) &dev_global, dim_globals);
    
    printf("Init coords\n");
    strcpy(shot.name, "17178_1097");
    init_taiga_props("particles", "1000", &beam, &shot, &run);
    
    printf("Init coords\n");
    init_coords(beam, shot, run, host_global, shared_global);
    set_particle_number(host_global, &run);
    set_particle_number(shared_global, &run);
    cudaMemcpy(dev_global, shared_global, dim_globals, cudaMemcpyHostToDevice);
        
    double *h_tmp, *d_tmp;
    size_t dim_tmp = sizeof(double)*LENGTH_TMP;
    h_tmp = (double *) malloc(dim_tmp);
    init_tmp(h_tmp);
    cudaMalloc((void **) &(d_tmp), dim_tmp);
    cudaMemcpy(d_tmp, h_tmp, dim_tmp, cudaMemcpyHostToDevice);
    
    test_init_coords_cuda <<< 1, 1 >>> (dev_global, d_tmp);
    
    cudaMemcpy(h_tmp, d_tmp, dim_tmp, cudaMemcpyDeviceToHost);
    
    print_tmp(h_tmp);
    
    printf("Number of ions: %d\n", host_global->particle_number);
    
    printf("1st:  [%lf, %lf, %lf]\n", host_global->rad[0], host_global->z[0], host_global->tor[0]);
    printf("2nd:  [%lf, %lf, %lf]\n", host_global->rad[1], host_global->z[1], host_global->tor[1]);
    printf("Last: [%lf, %lf, %lf]\n", host_global->rad[host_global->particle_number-1], host_global->z[host_global->particle_number-1], host_global->tor[host_global->particle_number-1]);
    
}

int main(){
    test_init_grid();
    test_init_coords();
}
