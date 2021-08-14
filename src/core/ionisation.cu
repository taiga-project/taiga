#include "ionisation.cuh"

__device__ void no_ionisation_loss(double *X, TaigaCommons *c){
    ;
}

__device__ void calculate_ionisation_loss(double *X, TaigaCommons *c){
    //get_temperature()
    //get_density()
    //X[BEAMLET_INTENSITY_ID] -= timestep * calculate_ioniation_rate(temperature, density, ionisation_energy);
}

//NRL Plasma Formulary 2016, p55
__device__ double calculate_ioniation_rate(double temperature, double density, double ionisation_energy){
    double relative_temperature = temperature / ionisation_energy;
    return 1e-11 * density * sqrt(relative_temperature) /
           (ionisation_energy * sqrt(ionisation_energy) * 6.0 + relative_temperature) *
           exp( -1.0/relative_temperature);
}
