#ifndef LORENTZ_CUH
#define LORENTZ_CUH

__device__ void get_acceleration_from_lorentz_force_with_electric_field(double *a, double *v,
                                                                        double *B, double *E,
                                                                        double eperm);

__device__ void get_acceleration_from_lorentz_force_without_electric_field(double *a, double *v,
                                                                           double *B, double *E,
                                                                           double eperm);

#endif //LORENTZ_CUH
