#include "helper.cuh"

#define LENGTH_TMP 10

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
