#include "core/init/init_beamlet.cuh"

template<typename Solver, typename MagneticField,
         typename SecondaryIonisation,
         typename Lorentz,
         typename DetectInterp,
         typename Perturb>
__global__ void taiga(TaigaGlobals *g, TaigaCommons *c, double *service_var) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (g->detcellid[idx] == -1) {
        double X[X_SIZE];
        int detcellid = g->detcellid[idx];

        // Load state
        X[0] = g->rad[idx];
        X[1] = g->z[idx];
        X[2] = g->tor[idx];
        X[3] = g->vrad[idx];
        X[4] = g->vz[idx];
        X[5] = g->vtor[idx];
        X[BEAMLET_INTENSITY_ID] = g->intensity[idx];
        X[TIME_OF_FLIGHT_ID] = g->time_of_flight[idx];

        detcellid = calculate_trajectory
                <Solver, MagneticField, SecondaryIonisation, Lorentz, Perturb>
        (c, X, detcellid);

        // Store back
        g->detcellid[idx]      = detcellid;
        g->rad[idx]            = X[0];
        g->z[idx]              = X[1];
        g->tor[idx]            = X[2];
        g->vrad[idx]           = X[3];
        g->vz[idx]             = X[4];
        g->vtor[idx]           = X[5];
        g->intensity[idx]      = X[BEAMLET_INTENSITY_ID];
        g->time_of_flight[idx] = X[TIME_OF_FLIGHT_ID];
    }

    // Example service var write
    if (threadIdx.x == 0 && blockIdx.x == 0) service_var[0] = 42.24;
}
