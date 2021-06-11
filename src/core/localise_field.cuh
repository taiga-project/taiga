#ifndef LOCALISE_FIELD_CUH
#define LOCALISE_FIELD_CUH

__device__ void get_coefficients_with_splines(
                                 TaigaCommons *c,
                                 const int *local_spline_indices,
                                 double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                 double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                 double *local_polflux, int zgrid_length);

__device__ void copy_local_field(TaigaCommons *c,
                                 double position_rad, double position_z,
                                 int *local_spline_indices,
                                 double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                 double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                 double *local_polflux);

__device__ double calculate_local_field_with_splines(TaigaCommons *c, const int *local_spline_indices,
                                                     const double *local_spline, double dr, double dz);

__device__ double calculate_local_field_with_bsplines(TaigaCommons *c, const int *local_spline_indices,
                                                      const double *local_spline, double dr, double dz);

#endif //LOCALISE_FIELD_CUH
