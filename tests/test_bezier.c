#include "test_bezier.h"
#include "taiga_test.h"
#include "utils/taiga_constants.h"
#define __device__ ;
#include "core/maths/maths.cu"

double get_bezier(double X_prev[6], double X[6], double detector[4], double timestep, int index){
    interpolate_bezier(X_prev, X, detector, timestep);
    return X[index];
}

int test_bezier() {
    TAIGA_INIT_TEST("BEZIER");

    double X_prev[6] = {-1, 1, 0, sqrt(2), sqrt(2), 0};
    double X[6] = {2, 0, 0, 0, -1, 0};
    double D[4] = {1, 0, 0, -1};
    TAIGA_ASSERT_ALMOST_EQ(1, get_bezier(X_prev, X, D, 1, 0), "x @ x=1 vs x^2+y^2=4");

    double X2_prev[6] = {0, sqrt(3), 0, sqrt(3), 1, 0};
    double X2[6] = {2, sqrt(3), 0, sqrt(3), -1, 0};
    TAIGA_ASSERT_ALMOST_EQ(1, get_bezier(X2_prev, X2, D, 1, 0), "x @ x=1 vs (x-1)^2+y^2=4");

    double X3_prev[6] = {1, 2, 0, 1, 1, 0};
    double X3[6] = {2, 1, 0, 1, -1, 0};
    double D3[4] = {1.0, 1.0, 0, -1.0};
    TAIGA_ASSERT_ALMOST_EQ(2, get_bezier(X3_prev, X3, D3, 1, 0), "x @ x=1 vs (x-1)^2+y^2=4");
    TAIGA_ASSERT_ALMOST_EQ(1, get_bezier(X3_prev, X3, D3, 1, 1), "y @ x=1 vs (x-1)^2+y^2=4");


    return TAIGA_ASSERT_SUMMARY();
}