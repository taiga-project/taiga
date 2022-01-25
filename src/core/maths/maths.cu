#include "core/maths/maths.cuh"

__device__ double cross(double *u, double *v, int index) {
    switch(index){
        case 0:
            return u[1]*v[2] - u[2]*v[1];
        case 1:
            return u[2]*v[0] - u[0]*v[2];
        case 2:
            return u[0]*v[1] - u[1]*v[0];
    }
}

__device__ double interpolate(double y1, double y2, double x, double x1, double x2){
    double x2_minus_x1 = x2-x1;
    return (x2-x) / x2_minus_x1 * y1 + (x - x1) / x2_minus_x1 * y2;
}

__device__ double interpolate_from_vector(double *x_vector, double *y_vector, long length, double x_value){
    long i=0;
    for (i=0; (i<length-1) && (x_vector[i]>x_value); ++i);
    return y_vector[i+1] - (y_vector[i+1]-y_vector[i])*(x_value-x_vector[i])/(x_vector[i+1]-x_vector[i]);
}