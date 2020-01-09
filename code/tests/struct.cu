#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "../prop.h"

#define LENGTH_DETCELLID 3
#define LENGTH_TMP 10
#define DETCELLID_INDEX 0

__global__ void test_struct_copy(device_global* g, double *tmp){
    tmp[0] = 3.14159265358979324;
    tmp[1] = g->detcellid[DETCELLID_INDEX];
    g->particle_number = 20;    
    g->detcellid[DETCELLID_INDEX] = 24;
    tmp[2] = g->detcellid[DETCELLID_INDEX];
}

void print_tmp(double *h_tmp){
    for (int i=0; i<LENGTH_TMP; ++i){
        printf("TMP %d: %lf\n", i, h_tmp[i]);
    }
}

int main(){
    device_global *d_global, *h_global, *s_global;
    size_t dim_global = sizeof(device_global);
    size_t dim_detcellid = sizeof(int)*LENGTH_DETCELLID;
    
    h_global = (device_global*)malloc(dim_global);   
    s_global = (device_global*)malloc(dim_global);    
    cudaMalloc((void **) &d_global, dim_global);
    
    h_global->detcellid = (int*)malloc(dim_detcellid);
        
    int *s_global__detcellid;
    cudaMalloc((void **) &(s_global__detcellid), dim_detcellid);
    
    h_global->particle_number = 10;
    h_global->detcellid[DETCELLID_INDEX] = 42;
    
    memcpy(s_global, h_global, dim_global);
    cudaMemcpy(s_global__detcellid, h_global->detcellid, dim_detcellid, cudaMemcpyHostToDevice);
    s_global->detcellid = s_global__detcellid;
    cudaMemcpy(d_global, s_global, dim_global, cudaMemcpyHostToDevice);
   
    double *h_tmp, *d_tmp;
    size_t dim_tmp = sizeof(double)*LENGTH_TMP;
    h_tmp = (double *) malloc(dim_tmp);
    cudaMalloc((void **) &(d_tmp), dim_tmp);
    cudaMemcpy(d_tmp, h_tmp, dim_tmp, cudaMemcpyHostToDevice);
    
    test_struct_copy <<< 1, 1 >>> (d_global, d_tmp);
    
    device_global *h2_global, *s2_global;
    int *s2_global__detcellid;
    s2_global__detcellid = (int*)malloc(dim_detcellid);    
     
    h2_global = (device_global*)malloc(dim_global);
    s2_global = (device_global*)malloc(dim_global);
    h2_global->detcellid = (int*)malloc(dim_detcellid);
    
    cudaMemcpy(s2_global, d_global, dim_global, cudaMemcpyDeviceToHost);
    cudaMemcpy(s2_global__detcellid, s2_global->detcellid, dim_detcellid, cudaMemcpyDeviceToHost);
    memcpy(h2_global, s2_global, dim_global);
    h2_global->detcellid = s2_global__detcellid;
    
    cudaMemcpy(h_tmp, d_tmp, dim_tmp, cudaMemcpyDeviceToHost);
    
    print_tmp(h_tmp);
    
    printf("TEST: Particle number (10): %d\n", h_global->particle_number);
    printf("TEST: Particle number (20): %d\n", h2_global->particle_number);
    printf("TEST: Detcellid[%d] (42): %d\n", DETCELLID_INDEX, h_global->detcellid[DETCELLID_INDEX]);
    printf("TEST: Detcellid[%d] (24): %d\n", DETCELLID_INDEX, h2_global->detcellid[DETCELLID_INDEX]);
    
}
