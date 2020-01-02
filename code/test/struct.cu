#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "../prop.h"

__global__ void test_struct_copy(device_global* g){
    g->particle_number = 20;
}

int main(){  
    device_global *d_global, *h_global, *h2_global;
    size_t dim_global = sizeof(device_global);
    
    h_global = (device_global*)malloc(dim_global);    
    h2_global = (device_global*)malloc(dim_global);
    cudaMalloc((void **) &d_global, dim_global);
    
    h_global->particle_number = 10;
    
    cudaMemcpy(d_global, h_global, dim_global, cudaMemcpyHostToDevice);
    
    test_struct_copy <<< 1, 1 >>> (d_global);
    cudaMemcpy(h2_global, d_global, dim_global, cudaMemcpyDeviceToHost);
    
    printf("TEST: Particle number (10): %d\n", h_global->particle_number);
    printf("TEST: Particle number (20): %d\n", h2_global->particle_number);
    
}
