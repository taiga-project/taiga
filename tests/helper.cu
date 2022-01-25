#include <stdio.h>
#include "helper.cuh"
#include "utils/taiga_constants.h"
#include "taiga_test.h"

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

void test_tmp(int number_of_cases, double *ref_tmp, double *h_tmp){
    for (int i=0; i<number_of_cases; ++i){
        char t[10];
        sprintf(t, "test %d", i);
        TAIGA_ASSERT_EQ(ref_tmp[i], h_tmp[i], t);
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
}
