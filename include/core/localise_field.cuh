#ifndef LOCALISE_FIELD_CUH
#define LOCALISE_FIELD_CUH

__device__ double (*calculate_local_field)(TaigaCommons *c, const int *local_spline_indices,
                                           const double *local_spline, double dr, double dz);

__device__ double (*get_dr)(TaigaCommons *c, const int *local_spline_indices, double R);
__device__ double (*get_dz)(TaigaCommons *c, const int *local_spline_indices, double Z);

__device__ void get_coefficients_with_splines(
        TaigaCommons *c,
        int *local_spline_indices,
        double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
        double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
        double *local_psi_n, int rgrid_length, int zgrid_length);

__device__ void get_coefficients_with_bsplines(
        TaigaCommons *c,
        int *local_spline_indices,
        double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
        double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
        double *local_psi_n, int rgrid_length, int zgrid_length);

__device__ double get_local_field(double *X, double *local_bfield, double *local_efield,
                                  TaigaCommons *c, bool is_electric_field_on,
                                  int *local_spline_indices,
                                  double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                  double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                  double *local_spline_psi_n);

__device__ void copy_local_field_coefficients(TaigaCommons *c,
                                              double R, double *X,
                                              int *local_spline_indices,
                                              double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                              double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                              double *local_spline_psi_n);

__device__ double get_local_field_from_coefficients(double *local_bfield, double *local_efield, bool is_electric_field_on,
                                                    TaigaCommons *c,
                                                    double R, double *X,
                                                    int *local_spline_indices,
                                                    double *local_spline_brad, double *local_spline_bz, double *local_spline_btor,
                                                    double *local_spline_erad, double *local_spline_ez, double *local_spline_etor,
                                                    double *local_spline_psi_n);

__device__ double calculate_local_field_with_splines(TaigaCommons *c, const int *local_spline_indices,
                                                     const double *local_spline, double dr, double dz);

__device__ double calculate_local_field_with_bsplines(TaigaCommons *c, const int *local_spline_indices,
                                                      const double *local_spline, double dr, double dz);

__device__ double get_dr_with_splines(TaigaCommons *c, const int *local_spline_indices, double R);
__device__ double get_dz_with_splines(TaigaCommons *c, const int *local_spline_indices, double Z);
__device__ double get_dr_with_bsplines(TaigaCommons *c, const int *local_spline_indices, double R);
__device__ double get_dz_with_bsplines(TaigaCommons *c, const int *local_spline_indices, double Z);

#endif //LOCALISE_FIELD_CUH
