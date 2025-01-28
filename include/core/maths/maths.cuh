#ifndef MATHS_CUH
#define MATHS_CUH

__device__ double cross(double *u, double *v, int index);
__device__ double interpolate(double y1, double y2, double x, double x1, double x2);
__device__ double interpolate_from_vector(double *x_vector, double *y_vector, long length, double x_value);

__device__ double solve_quadratic(double a, double b, double c);
__device__ void interpolate_bezier(double X_prev[6], double X[6], double detector[4], double timestep);

#endif //MATHS_CUH
