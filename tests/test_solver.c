#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#include "taiga_test.h"

#include "test_solver.h"

#include "utils/taiga_constants.h"
#include "utils/prop.c"
#include "core/maths/maths.cu"
#include "core/solvers/rk4.cu"
#include "core/solvers/runge_kutta_nystrom.cu"
#include "core/solvers/yoshida.cu"
#include "core/solvers/verlet.cu"
#include "core/solvers/boris.cu"
#include "core/physics/lorentz.cu"

#define GRAD_B_FACTOR 0.01
#define E_OVER_B 1

#define NUMBER_OF_CYCLOTRON_PERIODS 1000
#define NUMBER_OF_CYCLOTRON_PERIODS_GRAD_B 10000
#define LARMOR_RADIUS 0.01

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

void generate_homogeneous_magnetic_field(double *X, double *local_bfield, double *local_efield,
                                         TaigaCommons *c, bool is_electric_field_on,
                                         int *local_spline_indices,
                                         double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                         double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                         double *local_psi_n) {
    local_bfield[0] = 0.0;
    local_bfield[1] = 0.0;
    local_bfield[2] = 1.0;
    local_efield[0] = 0.0;
    local_efield[1] = 0.0;
    local_efield[2] = 0.0;
}

void generate_homogeneous_electric_field(double *X, double *local_bfield, double *local_efield,
                                         TaigaCommons *c, bool is_electric_field_on,
                                         int *local_spline_indices,
                                         double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                         double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                         double *local_psi_n) {
    local_bfield[0] = 0.0;
    local_bfield[1] = 0.0;
    local_bfield[2] = 0.0;
    local_efield[0] = 0.0;
    local_efield[1] = 0.0;
    local_efield[2] = 1.0;
}

void generate_grad_B_field(double *X, double *local_bfield, double *local_efield,
                           TaigaCommons *c, bool is_electric_field_on,
                           int *local_spline_indices,
                           double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                           double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                           double *local_psi_n) {
    local_bfield[0] = 0.0;
    local_bfield[1] = 0.0;
    local_bfield[2] = 1.0 + GRAD_B_FACTOR * X[1];
    local_efield[0] = 0.0;
    local_efield[1] = 0.0;
    local_efield[2] = 0.0;
}

void generate_E_par_B_field(double *X, double *local_bfield, double *local_efield,
                                    TaigaCommons *c, bool is_electric_field_on,
                                    int *local_spline_indices,
                                    double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                    double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                    double *local_psi_n) {
    local_bfield[0] = 0.0;
    local_bfield[1] = 0.0;
    local_bfield[2] = 1.0;
    local_efield[0] = 0.0;
    local_efield[1] = 0.0;
    local_efield[2] = E_OVER_B;
}

void generate_inv_R_field(double *X, double *local_bfield, double *local_efield,
                          TaigaCommons *c, bool is_electric_field_on,
                          int *local_spline_indices,
                          double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                          double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                          double *local_psi_n) {
    local_bfield[0] = 0.0;
    local_bfield[1] = 0.0;
    local_bfield[2] = LARMOR_RADIUS / sqrt(X[0] * X[0] + X[1] * X[1]);
    local_efield[0] = 0.0;
    local_efield[1] = 0.0;
    local_efield[2] = 0.0;
}

void generate_R_field(double *X, double *local_bfield, double *local_efield,
                          TaigaCommons *c, bool is_electric_field_on,
                          int *local_spline_indices,
                          double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                          double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                          double *local_psi_n) {
    local_bfield[0] = 0.0;
    local_bfield[1] = 0.0;
    local_bfield[2] = sqrt(X[0] * X[0] + X[1] * X[1]) / LARMOR_RADIUS;
    local_efield[0] = 0.0;
    local_efield[1] = 0.0;
    local_efield[2] = 0.0;
}

double get_local_field(double *X, double *local_bfield, double *local_efield,
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

double run_field_with_solver(double timestep, int field_type, int return_type,
                             double (*solve_diffeq)(double *X, double eperm, double timestep,
                                                  TaigaCommons *c, bool is_electric_field_on,
                                                  int *local_spline_indices,
                                                  double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                                  double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                                  double *local_psi_n) ) {


    int maximum_extrema = 2 * NUMBER_OF_CYCLOTRON_PERIODS + 1;

    switch(field_type){
        case HOMOGENEOUS:
            generate_local_field = &generate_homogeneous_magnetic_field;
            break;
        case E_FIELD:
            generate_local_field = &generate_homogeneous_electric_field;
            break;
        case GRAD_B:
            generate_local_field = &generate_grad_B_field;
            maximum_extrema = 2 * NUMBER_OF_CYCLOTRON_PERIODS_GRAD_B + 1;
            break;
        case E_PAR_B:
            generate_local_field = &generate_E_par_B_field;
            break;
        case INV_R:
            generate_local_field = &generate_inv_R_field;
            break;
        case PROP_R:
            generate_local_field = &generate_inv_R_field;
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
    int is_electric_field_on = true;

    get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_with_electric_field;

    SolverTestExtrema t;
    t.index = 0;
    t.direction = 1;
    t.counter = 1;
    t.extrema = (double*)malloc(maximum_extrema * sizeof(double));

    X[4] = eperm * LARMOR_RADIUS;
    t.extrema[0] = X[t.index];

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

    switch (return_type) {
        case GET_POSITION:
            return t.extrema[maximum_extrema-1];
        case GET_SPEED:
            return get_speed(X);
        case GET_SPEED_TOROIDAL:
            return(X[5]);
    }

    return -99999999999;
}

void export_field_with_solver(double timestep, int field_type,
                              double (*solve_diffeq)(double *X, double eperm, double timestep,
                                                     TaigaCommons *c, bool is_electric_field_on,
                                                     int *local_spline_indices,
                                                     double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                                     double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                                     double *local_psi_n),
                              long number_of_periods, char* filename ) {

    switch(field_type){
        case HOMOGENEOUS:
            generate_local_field = &generate_homogeneous_magnetic_field;
            break;
        case E_FIELD:
            generate_local_field = &generate_homogeneous_electric_field;
            break;
        case GRAD_B:
            generate_local_field = &generate_grad_B_field;
            break;
        case E_PAR_B:
            generate_local_field = &generate_E_par_B_field;
            break;
        case INV_R:
            generate_local_field = &generate_inv_R_field;
            break;
        case PROP_R:
            generate_local_field = &generate_R_field;
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
    int is_electric_field_on = true;

    double steps_in_a_period = (2.0 * 3.14159265358979 / eperm / timestep)*100;
    long maximum_step = number_of_periods * (long)steps_in_a_period;

    get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_with_electric_field;

    X[0] = -LARMOR_RADIUS;
    X[4] = eperm * LARMOR_RADIUS;

    // Open file in write mode
    FILE *file_ptr;
    file_ptr = fopen(filename, "w");
    long i;
    for(i=0; i < maximum_step; ++i){
        memcpy(&X_prev, &X, 6*sizeof(double));
        if (i % (long)steps_in_a_period == 0){
            fprintf(file_ptr, "%.12le\t%.12le\t%.12le\t%.12le\t%.12le\t%.12le\n", X[0], X[1], X[2], X[3], X[4], X[5]);
        }
        solve_diffeq(X, eperm, timestep,
                     c, is_electric_field_on,
                     local_spline_indices,
                     local_spline_brad, local_spline_bz, local_spline_btor,
                     local_spline_erad, local_spline_ez, local_spline_etor,
                     local_psi_n);
 }

    // Close the file
    fclose(file_ptr);
}

int test_solver() {
    TAIGA_INIT_TEST("SOLVER");
    double timestep = 1e-9;
    double eperm = 4e7;

    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_field_with_solver(timestep, HOMOGENEOUS, GET_POSITION, solve_diffeq_by_rk4), 1e-5, "4th order linearised Runge--Kutta");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_field_with_solver(timestep, HOMOGENEOUS, GET_POSITION, solve_diffeq_by_rkn), 1e-5, "4th order Runge--Kutta--Nystrom");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_field_with_solver(timestep, HOMOGENEOUS, GET_POSITION, solve_diffeq_by_verlet), 1e-5, "velocity-Verlet--Boris");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(0.0, run_field_with_solver(timestep, HOMOGENEOUS, GET_POSITION, solve_diffeq_by_yoshida), 1e-5, "Yoshida--Boris");

    double reference_grad_B_drift = -GRAD_B_FACTOR * PI * LARMOR_RADIUS * LARMOR_RADIUS * NUMBER_OF_CYCLOTRON_PERIODS_GRAD_B;
    double reference_speed = eperm * LARMOR_RADIUS;
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_grad_B_drift, run_field_with_solver(timestep, GRAD_B, GET_POSITION, solve_diffeq_by_rk4), 3e-5, "4th order linearised Runge--Kutta (grad B)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_grad_B_drift, run_field_with_solver(timestep, GRAD_B, GET_POSITION, solve_diffeq_by_rkn), 3e-5, "4th order linearised Runge--Kutta--Nystrom (grad B)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_grad_B_drift, run_field_with_solver(timestep, GRAD_B, GET_POSITION, solve_diffeq_by_verlet), 3e-5, "velocity-Verlet--Boris (grad B)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_grad_B_drift, run_field_with_solver(timestep, GRAD_B, GET_POSITION, solve_diffeq_by_yoshida), 3e-5, "Yoshida--Boris (grad B)");

    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_speed, run_field_with_solver(timestep, GRAD_B, GET_SPEED, solve_diffeq_by_rk4), 1000, "4th order linearised Runge--Kutta (grad B, speed)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_speed, run_field_with_solver(timestep, GRAD_B, GET_SPEED, solve_diffeq_by_rkn), 1000, "4th order linearised Runge--Kutta--Nystrom (grad B, speed)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_speed, run_field_with_solver(timestep, GRAD_B, GET_SPEED, solve_diffeq_by_verlet), 1, "velocity-Verlet--Boris (grad B, speed)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_speed, run_field_with_solver(timestep, GRAD_B, GET_SPEED, solve_diffeq_by_yoshida), 1, "Yoshida--Boris (grad B, speed)");

    double reference_speed_in_electric_field =  eperm * E_OVER_B * NUMBER_OF_CYCLOTRON_PERIODS * 2.0 * PI * LARMOR_RADIUS / reference_speed;
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_speed_in_electric_field, run_field_with_solver(timestep, E_PAR_B, GET_SPEED_TOROIDAL, solve_diffeq_by_rk4), 1, "4th order linearised Runge--Kutta (E || B)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_speed_in_electric_field, run_field_with_solver(timestep, E_PAR_B, GET_SPEED_TOROIDAL, solve_diffeq_by_rkn), 1, "4th order linearised Runge--Kutta-Nystrom (E || B)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_speed_in_electric_field, run_field_with_solver(timestep, E_PAR_B, GET_SPEED_TOROIDAL, solve_diffeq_by_verlet), 1, "velocity-Verlet--Boris (E || B)");
    TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(reference_speed_in_electric_field, run_field_with_solver(timestep, E_PAR_B, GET_SPEED_TOROIDAL, solve_diffeq_by_yoshida), 1, "Yoshida--Boris (E || B)");

    return TAIGA_ASSERT_SUMMARY();
}

void export_field() {
    long number_of_periods = 3000;
    double timestep = 1e-9;

    export_field_with_solver(timestep, HOMOGENEOUS, solve_diffeq_by_rk4, number_of_periods, "export_homo_rk4.txt");
    export_field_with_solver(timestep, HOMOGENEOUS, solve_diffeq_by_rkn, number_of_periods, "export_homo_rkn.txt");
    export_field_with_solver(timestep, HOMOGENEOUS, solve_diffeq_by_verlet, number_of_periods, "export_homo_verlet.txt");
    export_field_with_solver(timestep, HOMOGENEOUS, solve_diffeq_by_yoshida, number_of_periods, "export_homo_yoshida.txt");

    export_field_with_solver(timestep, E_PAR_B, solve_diffeq_by_rk4, number_of_periods, "export_eparb_rk4.txt");
    export_field_with_solver(timestep, E_PAR_B, solve_diffeq_by_rkn, number_of_periods, "export_eparb_rkn.txt");
    export_field_with_solver(timestep, E_PAR_B, solve_diffeq_by_verlet, number_of_periods, "export_eparb_verlet.txt");
    export_field_with_solver(timestep, E_PAR_B, solve_diffeq_by_yoshida, number_of_periods, "export_eparb_yoshida.txt");

    export_field_with_solver(timestep, GRAD_B, solve_diffeq_by_rk4, number_of_periods, "export_gradb_rk4.txt");
    export_field_with_solver(timestep, GRAD_B, solve_diffeq_by_rkn, number_of_periods, "export_gradb_rkn.txt");
    export_field_with_solver(timestep, GRAD_B, solve_diffeq_by_verlet, number_of_periods, "export_gradb_verlet.txt");
    export_field_with_solver(timestep, GRAD_B, solve_diffeq_by_yoshida, number_of_periods, "export_gradb_yoshida.txt");

    export_field_with_solver(timestep, INV_R, solve_diffeq_by_rk4, number_of_periods, "export_invr_rk4.txt");
    export_field_with_solver(timestep, INV_R, solve_diffeq_by_rkn, number_of_periods, "export_invr_rkn.txt");
    export_field_with_solver(timestep, INV_R, solve_diffeq_by_verlet, number_of_periods, "export_invr_verlet.txt");
    export_field_with_solver(timestep, INV_R, solve_diffeq_by_yoshida, number_of_periods, "export_invr_yoshida.txt");

    export_field_with_solver(timestep, PROP_R, solve_diffeq_by_rk4, number_of_periods, "export_r_rk4.txt");
    export_field_with_solver(timestep, PROP_R, solve_diffeq_by_rkn, number_of_periods, "export_r_rkn.txt");
    export_field_with_solver(timestep, PROP_R, solve_diffeq_by_verlet, number_of_periods, "export_r_verlet.txt");
    export_field_with_solver(timestep, PROP_R, solve_diffeq_by_yoshida, number_of_periods, "export_r_yoshida.txt");

}