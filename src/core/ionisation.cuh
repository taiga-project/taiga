#ifndef IONISATION_CUH
#define IONISATION_CUH

__device__ void calculate_ionisation_loss(double psi_n, double *X, TaigaCommons *c,
                                          int *local_ts_index, double *local_ts_psi, double timestep);
__device__ void no_ionisation_loss(double psi_n, double *X, TaigaCommons *c,
                                   int *local_ts_index, double *local_ts_psi, double timestep) {;}
__device__ double calculate_ionisation_rate(double temperature, double density, double ionisation_energy);
__device__ void localise_profile(double psi_n, double *X, TaigaCommons *c,
                                 int *local_ts_index, double *local_ts_psi);

#endif //IONISATION_CUH
