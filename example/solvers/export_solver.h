#ifndef EXPORT_SOLVER_CUH
#define EXPORT_SOLVER_CUH

#include <stdbool.h>
#include "utils/taiga_constants.h"
#include "utils/prop.h"

#define __device__ ;

#define BR_FIELD 90
#define B_OVER_R_FIELD 91

void export_coordinate (FILE *f, double *X);
void generate_R_B_const_field(double *X, double *local_bfield, double *local_efield,
                              TaigaCommons *c, bool is_electric_field_on,
                              int *local_spline_indices,
                              double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                              double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                              double *local_psi_n);
void generate_B_over_R_field(double *X, double *local_bfield, double *local_efield,
                             TaigaCommons *c, bool is_electric_field_on,
                             int *local_spline_indices,
                             double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                             double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                             double *local_psi_n);
void run_field_with_solver_and_export(char* scenario_name, double timestep, int field_type, char* solver_name,
                                      long number_of_cyclotron_periods, long frequency_of_export,
                                      double (*solve_diffeq)(double *X, double eperm, double timestep,
                                                             TaigaCommons *c, bool is_electric_field_on,
                                                             int *local_spline_indices,
                                                             double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                                             double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                                             double *local_psi_n) );
void run_scenario(char* scenario_name,double timestep, int field_type,
                  long number_of_cyclotron_periods, long frequency_of_export);

#endif //EXPORT_SOLVER_CUH
