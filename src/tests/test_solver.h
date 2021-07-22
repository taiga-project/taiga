#ifndef TEST_SOLVER_CUH
#define TEST_SOLVER_CUH

#include <stdbool.h>
#include "utils/taiga_constants.h"
#include "utils/prop.h"

#define __device__ ;

#define HOMOGENEOUS 0
#define GRAD_B 1

typedef struct SolverTestExtremaTag {
    double *extrema;
    int index;
    int direction;
    int counter;
}SolverTestExtrema;

double test_interp1(double x, double x0, double x1, double y0, double y1);
void init_homogeneous_field(double *X, double *B);
int get_direction(SolverTestExtrema *t, double *X, double *X_prev);
void get_extrema(SolverTestExtrema *t, double *X, double *X_prev);
double get_speed(double *X);
void (*generate_local_field)(double *X, double *local_bfield, double *local_efield,
                     TaigaCommons *c, bool is_electric_field_on,
                     int *local_spline_indices,
                     double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                     double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                     double *local_psi_n);
void generate_homogeneous_field(double *X, double *local_bfield, double *local_efield,
                                TaigaCommons *c, bool is_electric_field_on,
                                int *local_spline_indices,
                                double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                double *local_psi_n);
double run_field_with_solver(double timestep, int field_type,
                             void (*solve_diffeq)(double *X, double eperm, double timestep,
                                                  TaigaCommons *c, bool is_electric_field_on,
                                                  int *local_spline_indices,
                                                  double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                                  double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                                  double *local_psi_n) );

void test_solver();

#endif //TEST_SOLVER_CUH
