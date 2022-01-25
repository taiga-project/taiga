#include "core/physics/ionisation.cuh"

__device__ void calculate_ionisation_loss(double psi_n, double *X, TaigaCommons *c,
                                          int *local_ts_index, double *local_ts_psi, double timestep) {
    double temperature, density, loss;
    localise_profile(psi_n, c, local_ts_index, local_ts_psi);
    temperature = interpolate(c->ts_temperature[local_ts_index[0]], c->ts_temperature[local_ts_index[0]+1],
                              psi_n, local_ts_psi[0], local_ts_psi[1]);
    density = interpolate(c->ts_density[local_ts_index[0]], c->ts_density[local_ts_index[0]+1],
                              psi_n, local_ts_psi[0], local_ts_psi[1]);
    loss = timestep * calculate_ionisation_rate(temperature, density, c->ionisation_energy);
    if (isnan(loss) == false) {
        X[BEAMLET_INTENSITY_ID] *= (1.0 - loss);
    }
}

//NRL Plasma Formulary 2016, p55
__device__ double calculate_ionisation_rate(double temperature, double density, double ionisation_energy) {

    double relative_temperature = temperature / ionisation_energy;
    return 1e-11 * density * sqrt(relative_temperature) /
           (ionisation_energy * sqrt(ionisation_energy) * (6.0 + relative_temperature)) *
           exp( -1.0/relative_temperature) +
           8.75e-39*pow(temperature,-4.5)* density* density;
}

__device__ void localise_profile(double psi_n, TaigaCommons *c, int *local_ts_index, double *local_ts_psi) {

    int ts_length = c->ts_length;

    if (local_ts_index[0] == SPLINE_INDEX_ERROR) {
        local_ts_index[0] = 0;
    }

    if ( (local_ts_psi[0] == UNDEFINED_FLOAT) || (psi_n < local_ts_psi[0]) || (psi_n > local_ts_psi[1]) ) {
        while ((local_ts_index[0] > 0) && (psi_n < c->ts_psi[local_ts_index[0]])) {
            --local_ts_index[0];
        }
        while ((local_ts_index[0] < ts_length - 1) && (psi_n > c->ts_psi[local_ts_index[0] + 1])) {
            ++local_ts_index[0];
        }
        local_ts_psi[0] = c->ts_psi[local_ts_index[0]];
        if (local_ts_index[1] >= ts_length-1) {
            local_ts_psi[1] = INFINITY_INTEGER;
        }else {
            local_ts_psi[1] = c->ts_psi[local_ts_index[0] + 1];
        }
    }
}