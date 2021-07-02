#include <stdio.h>
#include <stdlib.h>

#define __device__ ;
#include "taiga_test.h"

#include "utils/taiga_constants.h"
#include "core/maths.cu"
#include "core/rk4.cu"
#include "core/yoshida.cu"
#include "core/verlet.cu"
#include "core/lorentz.cu"

struct SolverTestExtrema {
    double *extrema;
    int index;
    int direction;
    int counter;
};

inline double test_interp1(double x, double x0, double x1, double y0, double y1) {
    return y0+(x-x0)*(y1-y0)/(x1-x0);
}

void init_homogeneous_field(double *X, double *B) {
    X[4] = 4e5;
    B[2] = 1.0;
}

inline int get_direction(SolverTestExtrema *t, double *X, double *X_prev) {
    return ((X[t->index]-X_prev[t->index])>0) ? 1 : -1;
}

void get_extrema(SolverTestExtrema *t, double *X, double *X_prev) {
    int direction = get_direction(t, X, X_prev);
    if (direction != t->direction) {
        t->direction = direction;
        t->extrema[t->counter] = test_interp1(0, X_prev[3+t->index], X[3+t->index], X_prev[t->index], X[t->index]);
        t->counter++;
    }
}

double get_speed(double *X){
    return sqrt(X[3]*X[3]+X[4]*X[4]+X[5]*X[5]);
}

double run_homogeneous_field_with_solver(double timestep,
                                         void (*solve_diffeq)(double *X, double *B, double *E, double *E_prev,
                                                 double eperm, double timestep)){
    double X[6] = {0};
    double X_prev[6] = {0};
    double a[3] = {0};
    double B[3] = {0};
    double E[3] = {0};
    double E_prev[3] = {0};
    double eperm = 4e7;

    get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_without_electric_field;

    int maximum_extrema = 20;

    SolverTestExtrema t;
    t.index = 0;
    t.direction = -1;
    t.counter = 1;
    t.extrema = (double*)malloc(maximum_extrema * sizeof(double));

    init_homogeneous_field(X, B);
    t.extrema[0] = X[t.index];

    double v0 = get_speed(X);

    while (t.counter < maximum_extrema) {
        memcpy(&X_prev, &X, 6*sizeof(double));
        solve_diffeq(X, B, E, E_prev, eperm, timestep);
        get_extrema(&t, X, X_prev);
    }
    return t.extrema[maximum_extrema-1];
}

void test_solver() {
    TAIGA_INIT_TEST();
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_homogeneous_field_with_solver(1e-9, solve_diffeq_by_rk4), 1e-5, "4th order linearised Runge--Kutta");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_homogeneous_field_with_solver(1e-9, solve_diffeq_by_verlet), 1e-5, "velocity-Verlet based Boris-SDC (BGSDC)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_homogeneous_field_with_solver(1e-9, solve_diffeq_by_yoshida), 1e-5, "Yoshida based Boris-SDC");
    TAIGA_ASSERT_SUMMARY();
}

int main(){
    test_solver();
}