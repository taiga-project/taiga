#ifndef IONISATION_CUH
#define IONISATION_CUH

__device__ void calculate_ionisation_loss(double *X, TaigaCommons *c);
__device__ void no_ionisation_loss(double *X, TaigaCommons *c);
__device__ double calculate_ioniation_rate(double temperature, double density, double ionisation_energy);

#endif //IONISATION_CUH
