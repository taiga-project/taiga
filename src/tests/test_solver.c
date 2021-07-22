#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "taiga_test.h"

#include "test_solver.h"

#include "utils/taiga_constants.h"
#include "utils/prop.c"
#include "core/maths.cu"
#include "core/rk4.cu"
#include "core/runge_kutta_nystrom.cu"
#include "core/yoshida.cu"
#include "core/verlet.cu"
#include "core/lorentz.cu"

double test_interp1(double x, double x0, double x1, double y0, double y1) {
    return y0+(x-x0)*(y1-y0)/(x1-x0);
}

int get_direction(SolverTestExtrema *t, double *X, double *X_prev) {
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

double get_speed(double *X) {
    return sqrt(X[3]*X[3]+X[4]*X[4]+X[5]*X[5]);
}

void generate_homogeneous_field(double *X, double *local_bfield, double *local_efield,
                                TaigaCommons *c, bool is_electric_field_on,
                                int *local_spline_indices,
                                double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                double *local_psi_n){
    local_bfield[0] = 0.0;
    local_bfield[1] = 0.0;
    local_bfield[2] = 1.0;
    local_efield[0] = 0.0;
    local_efield[1] = 0.0;
    local_efield[2] = 0.0;
}

void get_local_field(double *X, double *local_bfield, double *local_efield,
                           TaigaCommons *c, bool is_electric_field_on,
                           int *local_spline_indices,
                           double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                           double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                           double *local_psi_n){
    generate_local_field(X, local_bfield, local_efield,
                         c, is_electric_field_on,
                         local_spline_indices,
                         local_spline_brad, local_spline_bz, local_spline_btor,
                         local_spline_erad, local_spline_ez, local_spline_etor,
                         local_psi_n);
}

double run_field_with_solver(double timestep, int field_type,
                             void (*solve_diffeq)(double *X, double eperm, double timestep,
                                                  TaigaCommons *c, bool is_electric_field_on,
                                                  int *local_spline_indices,
                                                  double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                                  double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                                  double *local_psi_n) ) {

    switch(field_type){
        case HOMOGENEOUS:
            generate_local_field = &generate_homogeneous_field;
            break;
        default:
            printf("Error: Illegal field_type\n");
    }
    double X[6] = {0};
    double X_prev[6] = {0};
    int local_spline_indices[2];
    local_spline_indices[0] = SPLINE_INDEX_ERROR;
    local_spline_indices[1] = SPLINE_INDEX_ERROR;
    double local_spline_brad[16];
    double local_spline_bz[16];
    double local_spline_btor[16];
    double local_spline_erad[16];
    double local_spline_ez[16];
    double local_spline_etor[16];
    double local_psi_n[16];
    TaigaCommons *c;
    double eperm = 4e7;
    int is_electric_field_on = false;

    get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_without_electric_field;

    int maximum_extrema = 20;

    SolverTestExtrema t;
    t.index = 0;
    t.direction = -1;
    t.counter = 1;
    t.extrema = (double*)malloc(maximum_extrema * sizeof(double));

    X[4] = 4e5;
    t.extrema[0] = X[t.index];

    double v0 = get_speed(X);

    while (t.counter < maximum_extrema) {
        memcpy(&X_prev, &X, 6*sizeof(double));
        solve_diffeq(X, eperm, timestep,
                     c, is_electric_field_on,
                     local_spline_indices,
                     local_spline_brad, local_spline_bz, local_spline_btor,
                     local_spline_erad, local_spline_ez, local_spline_etor,
                     local_psi_n);
        get_extrema(&t, X, X_prev);
    }
    return t.extrema[maximum_extrema-1];
}

void test_solver() {
    TAIGA_INIT_TEST("SOLVER");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_field_with_solver(1e-9, HOMOGENEOUS, solve_diffeq_by_rk4), 1e-5, "4th order linearised Runge--Kutta");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_field_with_solver(1e-9, HOMOGENEOUS, solve_diffeq_by_rkn), 1e-5, "4th order Runge--Kutta--Nystrom");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_field_with_solver(1e-9, HOMOGENEOUS, solve_diffeq_by_verlet), 1e-5, "velocity-Verlet based Boris-SDC (BGSDC)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_field_with_solver(1e-9, HOMOGENEOUS, solve_diffeq_by_yoshida), 1e-5, "Yoshida based Boris-SDC");
    TAIGA_ASSERT_SUMMARY();
}