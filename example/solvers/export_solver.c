#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

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

#define GRAD_B_FACTOR 0.0001
#define E_OVER_B 1
#define LARMOR_RADIUS 0.01

#define FOLDER "example/solvers/data"

void export_coordinate (FILE *f, double *X) {
    fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", X[0], X[1], X[2], X[3], X[4], X[5]);
}

void run_field_with_solver_and_export(double timestep, int field_type, char* file_name, long number_of_cyclotron_periods,
                                      double (*solve_diffeq)(double *X, double eperm, double timestep,
                                                             TaigaCommons *c, bool is_electric_field_on,
                                                             int *local_spline_indices,
                                                             double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                                             double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                                             double *local_psi_n) ) {

    int maximum_extrema = 2 * number_of_cyclotron_periods + 1;

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

    get_acceleration_from_lorentz_force = &get_acceleration_from_lorentz_force_with_electric_field;
    X[4] = eperm * LARMOR_RADIUS;

    file = fopen(concat(FOLDER, "/", file_name, ".dat", NULL) ,"w");
    for (int i = 0; i < maximum_extrema; ++i) {
        solve_diffeq(X, eperm, timestep,
                     c, is_electric_field_on,
                     local_spline_indices,
                     local_spline_brad, local_spline_bz, local_spline_btor,
                     local_spline_erad, local_spline_ez, local_spline_etor,
                     local_psi_n);
        export_coordinate(file, X);
    }
    fclose(file);
}

int main() {
    double timestep = 1e-12;
    long number_of_cyclotron_periods = 10000;
    run_field_with_solver_and_export(timestep, HOMOGENEOUS, "rk4", number_of_cyclotron_periods, solve_diffeq_by_rk4);
    run_field_with_solver_and_export(timestep, HOMOGENEOUS, "rkn", number_of_cyclotron_periods, solve_diffeq_by_rkn);
    run_field_with_solver_and_export(timestep, HOMOGENEOUS, "verlet", number_of_cyclotron_periods,solve_diffeq_by_verlet);
    run_field_with_solver_and_export(timestep, HOMOGENEOUS, "yoshida", number_of_cyclotron_periods, solve_diffeq_by_yoshida);
    return 0;
}
