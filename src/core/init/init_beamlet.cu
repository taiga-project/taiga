#include "core/init/init_beamlet.cuh"
/*
__device__ void init_beamlet_intensity(double X[X_SIZE], int coordinate_id) {
    X[coordinate_id] = 1.0;
}

__device__ void init_time_of_flight(double X[X_SIZE], int coordinate_id,  TaigaCommons *c) {
    X[coordinate_id] = (X[0] - c->beam.deflection_radial_coordinate) / X[3];
}
*/