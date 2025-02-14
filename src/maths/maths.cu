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

// Bezier

__device__ double solve_quadratic(double a, double b, double c) {
    double D = b * b - 4.0 * a * c;
    if (D >= 0) {
        double sqrtD__2a = sqrt(D) / 2.0 / a;
        double mb__2a = -b / 2.0 / a;
        double t = mb__2a + sqrtD__2a;
        if ( t>= 0 && t <= 1) {
            return t;
        }
        t = mb__2a - sqrtD__2a;
        if ( t>= 0 && t <= 1) {
            return t;
        }
    }
    return OUT_OF_RANGE;
}

__device__ void interpolate_bezier(double X_prev[6], double X[6], double D[4], double timestep) {
    double P[3];

    for (int i = 0 ; i < 3; ++i) {
        P[i] = X_prev[i] + (1.0 / 3.0) * X_prev[i+3] * timestep;
    }

    double a = D[0] * (X_prev[0] - 2.0 * P[0] + X[0]) +
               D[1] * (X_prev[1] - 2.0 * P[1] + X[1]) +
               D[2] * (X_prev[2] - 2.0 * P[2] + X[2]);
    double b = 2.0 * D[0] * (P[0] - X_prev[0]) +
               2.0 * D[1] * (P[1] - X_prev[1]) +
               2.0 * D[2] * (P[2] - X_prev[2]);
    double c = D[0] * X_prev[0] + D[1] * X_prev[1] + D[2] * X_prev[2] + D[3];

    double t = solve_quadratic(a, b, c);

    if (t == OUT_OF_RANGE) {
        return;
    }

    for (int i = 0 ; i < 3; ++i) {
        X[i] = (1.0 - t) * (1.0 - t) * X_prev[i] +
                2.0 * (1.0 - t) * t * P[i] +
                t * t * X[i];
    }
}