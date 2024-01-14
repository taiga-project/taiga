#include "utils/free.cuh"

void free_taiga(TaigaGlobals *host_global, TaigaGlobals *shared_global, TaigaGlobals *device_global,
                TaigaCommons *host_common,  TaigaCommons *shared_common,  TaigaCommons *device_common,
                DetectorProp *shared_detector, DetectorProp *device_detector,
                ShotProp *shot, BeamProp *beam, RunProp *run,
                double *host_service_array, double *device_service_array) {
    if (run->debug == 1) printf("Free global.\n");
    free_global(host_global, HOST);
    free_global(shared_global, SHARED);
    free_global(device_global, DEVICE);
    if (run->debug == 1) printf("Free common.\n");
    free_common(host_common, HOST);
    free_common(shared_common, SHARED);
    free_common(device_common, DEVICE);
    if (run->debug == 1) printf("Free detector.\n");
    free_detector(shared_detector, SHARED);
    free_detector(device_detector, DEVICE);
    if (run->debug == 1) printf("Free service array.\n");
    free(host_service_array);
    CHECK_ERROR(cudaFree(device_service_array));
}

void free_detector(DetectorProp *detector, int where){
    FREE(where, detector->counter);
    FREE(where, detector->xgrid);
    FREE(where, detector->ygrid);
    FREE_MAIN(where, detector);
}

void free_beam(BeamProfile *beam_profile, int where){
    FREE(where, beam_profile->radial_grid);
    FREE(where, beam_profile->radial_profile);
    FREE(where, beam_profile->cross_grid);
    FREE(where, beam_profile->cross_profile);
    FREE_MAIN(where, beam_profile);
}

void free_global(TaigaGlobals *global, int where){
    FREE(where, global->rad);
    FREE(where, global->z);
    FREE(where, global->tor);
    FREE(where, global->vrad);
    FREE(where, global->vz);
    FREE(where, global->vtor);
    FREE(where, global->detcellid);
    FREE(where, global->intensity);
    FREE(where, global->time_of_flight);
    FREE_MAIN(where, global);
}

void free_common(TaigaCommons *common, int where){
    FREE(where, common->grid_size);
    FREE(where, common->spline_rgrid);
    FREE(where, common->spline_zgrid);
    FREE_SHARED(where, common->brad);
    FREE_SHARED(where, common->bz);
    FREE_SHARED(where, common->btor);
    FREE_SHARED_IF(where, common->erad);
    FREE_SHARED_IF(where, common->ez);
    FREE_SHARED_IF(where, common->etor);
    FREE_SHARED_IF(where, common->psi_n);
    FREE(where, common->detector_geometry);
    FREE_SHARED(where, common->ts_psi);
    FREE_SHARED(where, common->ts_density);
    FREE_SHARED(where, common->ts_temperature);
    FREE_MAIN(where, common);
}
