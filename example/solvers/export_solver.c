#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include "export_solver.h"
#include "test_solver.h"

#include "utils/taiga_constants.h"
#include "utils/prop.h"
#include "core/maths/maths.cuh"
#include "core/solvers/solvers.cuh"
#include "core/solvers/rk4.cuh"
#include "core/solvers/runge_kutta_nystrom.cuh"
#include "core/solvers/yoshida.cuh"
#include "core/solvers/verlet.cuh"
#include "core/solvers/boris.cuh"
#include "core/physics/lorentz.cuh"

#include "utils/basic_functions.c"

#define FOLDER "example"

void export_coordinate (FILE *f, double *X, long N_half) {
    fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%ld\n", X[0], X[1], X[2], X[3], X[4], X[5], N_half);
}

void generate_R_B_const_field(double *X, double *local_bfield, double *local_efield,
                              TaigaCommons *c, bool is_electric_field_on,
                              int *local_spline_indices,
                              double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                              double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                              double *local_psi_n) {
    local_bfield[0] = 0.0;
    local_bfield[1] = 0.0;
    local_bfield[2] = 1.0 / (X[0] + 1.0);
    local_efield[0] = 0.0;
    local_efield[1] = 0.0;
    local_efield[2] = 0.0;
}

void generate_B_over_R_field(double *X, double *local_bfield, double *local_efield,
                             TaigaCommons *c, bool is_electric_field_on,
                             int *local_spline_indices,
                             double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                             double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                             double *local_psi_n) {
    local_bfield[0] = 0.0;
    local_bfield[1] = 0.0;
    local_bfield[2] = sqrt (X[0] * X[0] + X[1] * X[1]);
    local_efield[0] = 0.0;
    local_efield[1] = 0.0;
    local_efield[2] = 0.0;
}

void run_field_with_solver_and_export(char* scenario_name, double timestep, int field_type, char* solver_name,
                                      long number_of_steps, long frequency_of_export,
                                      double (*solve_diffeq)(double *X, double eperm, double timestep,
                                                             TaigaCommons *c, bool is_electric_field_on,
                                                             int *local_spline_indices,
                                                             double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                                             double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                                             double *local_psi_n) ) {
    switch(field_type){
        case HOMOGENEOUS:
            generate_local_field = &generate_homogeneous_magnetic_field;
            break;
        case GRAD_B:
            generate_local_field = &generate_grad_B_field;
            break;
        case E_PAR_B:
            generate_local_field = &generate_E_par_B_field;
            break;
        case BR_FIELD:
            generate_local_field = &generate_R_B_const_field;
            break;
        case B_OVER_R_FIELD:
            generate_local_field = &generate_B_over_R_field;
            break;
        default:
            printf("Error: Illegal field_type\n");
    }
    FILE *file;
    double X[6] = {0};
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
    long number_of_half_cyclotron_periods = 0;
    double v_previous;

    get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_with_electric_field;
    X[4] = eperm * LARMOR_RADIUS;

    if (field_type == B_OVER_R_FIELD)   X[0] = 1.5;

    mkdir(FOLDER, S_IRWXU | S_IRWXG | S_IRWXO);
    const char *path = concat(FOLDER, "/", solver_name, "_", scenario_name, ".dat", NULL);
    printf("export to: %s\n", path);
    file = fopen(path ,"w");
    for (int i = 0; i < number_of_steps; ++i) {
        v_previous = X[4];
        solve_diffeq(X, eperm, timestep,
                     c, is_electric_field_on,
                     local_spline_indices,
                     local_spline_brad, local_spline_bz, local_spline_btor,
                     local_spline_erad, local_spline_ez, local_spline_etor,
                     local_psi_n);
        if ((X[4] * v_previous) <= 0)   ++number_of_half_cyclotron_periods;
        if (i % frequency_of_export == 0)   export_coordinate(file, X, number_of_half_cyclotron_periods);
    }
    fclose(file);
}


void run_scenario(char* scenario_name, double timestep, int field_type,
                  long number_of_steps, long frequency_of_export){
    run_field_with_solver_and_export(scenario_name, timestep, field_type, "rk4", number_of_steps, frequency_of_export, solve_diffeq_by_rk4);
    run_field_with_solver_and_export(scenario_name, timestep, field_type, "rkn", number_of_steps, frequency_of_export, solve_diffeq_by_rkn);
    run_field_with_solver_and_export(scenario_name, timestep, field_type, "verlet", number_of_steps, frequency_of_export, solve_diffeq_by_verlet);
    run_field_with_solver_and_export(scenario_name, timestep, field_type, "yoshida", number_of_steps, frequency_of_export, solve_diffeq_by_yoshida);
}

int main() {
    double timestep = 1e-9;
    long number_of_steps = 20000000;
    long frequency_of_export = 100000;

    run_scenario("default", timestep, HOMOGENEOUS, number_of_steps, frequency_of_export);
    run_scenario("gradb", timestep, GRAD_B, number_of_steps, frequency_of_export);
    run_scenario("eparb", timestep, E_PAR_B, number_of_steps, frequency_of_export);
    run_scenario("br", timestep, BR_FIELD, number_of_steps, frequency_of_export);
    run_scenario("b__r", timestep, B_OVER_R_FIELD, number_of_steps, frequency_of_export);
    run_scenario("start8", 1e-8, HOMOGENEOUS, 10000, 1);
    run_scenario("b__r_start8", 1e-8, B_OVER_R_FIELD, 10000, 1);
    run_scenario("start", timestep, HOMOGENEOUS, 100000, 10);
    run_scenario("b__r_start", timestep, B_OVER_R_FIELD, 100000, 10);
    run_scenario("gradb_start", timestep, GRAD_B, 100000, 10);
    run_scenario("eparb_start", timestep, E_PAR_B, 100000, 10);

    return 0;
}
