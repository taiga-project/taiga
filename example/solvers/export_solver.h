#ifndef EXPORT_SOLVER_CUH
#define EXPORT_SOLVER_CUH

#include <stdbool.h>
#include "utils/taiga_constants.h"
#include "utils/prop.h"

#define __device__ ;

void export_coordinate (FILE *f, double *X);
void run_field_with_solver_and_export(double timestep, int field_type, char* file_name, long number_of_cyclotron_periods,
                                      double (*solve_diffeq)(double *X, double eperm, double timestep,
                                                             TaigaCommons *c, bool is_electric_field_on,
                                                             int *local_spline_indices,
                                                             double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                                             double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                                             double *local_psi_n) );

#endif //EXPORT_SOLVER_CUH
