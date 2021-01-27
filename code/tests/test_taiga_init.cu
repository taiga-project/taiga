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
    init_grid(shot, run, host_common, shared_common, dev_common);
        
    double *h_tmp, *d_tmp;
    size_t dim_tmp = sizeof(double)*LENGTH_TMP;
    h_tmp = (double *) malloc(dim_tmp);
    init_tmp(h_tmp);
    cudaMalloc((void **) &(d_tmp), dim_tmp);
    cudaMemcpy(d_tmp, h_tmp, dim_tmp, cudaMemcpyHostToDevice);
    
    test_init_grid_cuda <<< 1, 1 >>> (dev_common, d_tmp);
    
    cudaMemcpy(h_tmp, d_tmp, dim_tmp, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    
    print_tmp(h_tmp);
    
    printf("R = %lf ... %lf\n", host_common->spline_rgrid[0], host_common->spline_rgrid[host_common->grid_size[0]-1]);
    printf("Z = %lf ... %lf\n", host_common->spline_zgrid[0], host_common->spline_zgrid[host_common->grid_size[1]-1]);
}

int main(){
    test_init_grid();
}
