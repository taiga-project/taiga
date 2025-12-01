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
    // Start and end positions
    double r1x = X_prev[0], r1y = X_prev[1], r1z = X_prev[2];
    double r2x = X[0],     r2y = X[1],     r2z = X[2];

    // Control point from derivative
    double dt_over3 = timestep / 3.0;
    double px = r1x + (1.0/3.0) * X_prev[3] * dt_over3;
    double py = r1y + (1.0/3.0) * X_prev[4] * dt_over3;
    double pz = r1z + (1.0/3.0) * X_prev[5] * dt_over3;

    // Quadratic Bézier projected onto plane normal D[0..2], offset D[3]
    double a = D[0]*(r1x - 2.0*px + r2x) +
               D[1]*(r1y - 2.0*py + r2y) +
               D[2]*(r1z - 2.0*pz + r2z);

    double b = 2.0 * (D[0]*(px - r1x) +
                      D[1]*(py - r1y) +
                      D[2]*(pz - r1z));

    double c = D[0]*r1x + D[1]*r1y + D[2]*r1z + D[3];

    // Solve quadratic for t
    double t = solve_quadratic(a, b, c);
    if (t == OUT_OF_RANGE) return;

    // Evaluate Bézier at t
    double omt = 1.0 - t;
    double omt2 = omt * omt;
    double t2   = t * t;

    X[0] = omt2 * r1x + 2.0 * omt * t * px + t2 * r2x;
    X[1] = omt2 * r1y + 2.0 * omt * t * py + t2 * r2y;
    X[2] = omt2 * r1z + 2.0 * omt * t * pz + t2 * r2z;
}


__device__ void interpolate_hermite(double X_prev[6], double X[6], double D[4], double timestep) {
    // Extract start and end positions
    double r1x = X_prev[0], r1y = X_prev[1], r1z = X_prev[2];
    double r2x = X[0],     r2y = X[1],     r2z = X[2];

    // Derivative at start scaled by timestep
    double v1x = X_prev[3] * timestep;
    double v1y = X_prev[4] * timestep;
    double v1z = X_prev[5] * timestep;

    // Quadratic Hermite coefficients
    double Ax = r2x - r1x - v1x;
    double Ay = r2y - r1y - v1y;
    double Az = r2z - r1z - v1z;

    double Bx = v1x;
    double By = v1y;
    double Bz = v1z;

    double Cx = r1x;
    double Cy = r1y;
    double Cz = r1z;

    // Project onto plane normal D[0..2], offset D[3]
    double a = D[0]*Ax + D[1]*Ay + D[2]*Az;
    double b = D[0]*Bx + D[1]*By + D[2]*Bz;
    double c = D[0]*Cx + D[1]*Cy + D[2]*Cz + D[3];

    // Solve quadratic for t
    double t = solve_quadratic(a, b, c);
    if (t == OUT_OF_RANGE) return;

    // Evaluate Hermite polynomial at t
    double t2 = t * t;
    X[0] = Ax*t2 + Bx*t + Cx;
    X[1] = Ay*t2 + By*t + Cy;
    X[2] = Az*t2 + Bz*t + Cz;
}
