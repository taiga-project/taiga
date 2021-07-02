#include <stdio.h>
#include <stdlib.h>
#define __device__ ;
#include "core/bspline.cu"
#include "taiga_test.h"

double test_bspline_scenario(double x0, double x1) {
    int i, j, index;
    int k = 3;
    int k_plus_1 = 4;

    int t0_length, t1_length;
    int i0, i1;

    double ref = 0.10051077;

    double t0[] = {-5., -5., -5., -5., -1., 1., 5., 5., 5., 5.};
    t0_length = 10;
    double t1[] = {-6., -6., -6., -6., -2., 0., 2., 6., 6., 6., 6.};
    t1_length = 11;

    double c[] = {0.9207306, 0.99199415, -2.58237284, 2.58237284, -0.99199415, -0.9207306,
                  -0.52069907, -0.56100062, 1.46040454, -1.46040454, 0.56100062, 0.52069907,
                  -1.72661944, -1.86025792, 4.84264905, -4.84264905, 1.86025792, 1.72661944,
                  1.85410551, 1.99761128, -5.2002092, 5.2002092, -1.99761128, -1.85410551,
                  -1.72661944, -1.86025792, 4.84264905, -4.84264905, 1.86025792, 1.72661944,
                  -0.52069907, -0.56100062, 1.46040454, -1.46040454, 0.56100062, 0.52069907,
                  0.9207306, 0.99199415, -2.58237284, 2.58237284, -0.99199415, -0.9207306};

    // copy_local_field()
    for (i0 = 0; (t0[i0 + 1] < x0) && (i0 < t0_length - 1); ++i0) { ; }
    for (i1 = 0; (t1[i1 + 1] < x1) && (i1 < t1_length - 1); ++i1) { ; }

    // get_coefficients_with_bsplines()
    if(i0<k_plus_1) i0=k;
    if(i1<k_plus_1) i1=k;

    double local_c[16];
    for (i = 0; i < k_plus_1; ++i) {
        for (j = 0; j < k_plus_1; ++j) {
            index = (i1 - k + j ) * (t0_length - k_plus_1)
                    + i0 - k + i;
            local_c[i*k_plus_1+j] = c[index];
        }
    }

    // calculate_local_field_with_bsplines()
    double B0[k_plus_1];
    double B1[k_plus_1];
    bspline(B0, x0, k, i0, t0);
    bspline(B1, x1, k, i1, t1);
    double  value = 0.0;
    for (i = 0; i < k_plus_1; ++i) {
        for (j = 0; j < k_plus_1; ++j) {
            value += local_c[i*k_plus_1+j] * B0[i] * B1[j];
        }
    }
    return  value;
}

void test_bspline() {
    TAIGA_INIT_TEST();
    TAIGA_ASSERT_ALMOST_EQ(0.10051077, test_bspline_scenario(0.1, 0.2), "bspline test 01");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, test_bspline_scenario(0.0, 0.0), 1e-7, "bspline test 02");
    TAIGA_ASSERT_ALMOST_EQ(0.95892427, test_bspline_scenario(-5.0, 0.0), "bspline test 03");
    TAIGA_ASSERT_ALMOST_EQ(0.95978649, test_bspline_scenario(-4.9, 0.0), "bspline test 04");
    TAIGA_ASSERT_SUMMARY();
}

int main(){
    test_bspline();
}

