#include "utils/free.cuh"




void free_taiga() {

}


void free_detector(DetectorProp *detector){
    CHECK_ERROR(cudaFree(detector->counter));
    CHECK_ERROR(cudaFree(detector->xgrid));
    CHECK_ERROR(cudaFree(detector->ygrid));
}

void free_beam(BeamProfile *beam_profile){
    CHECK_ERROR(cudaFree(beam_profile->radial_grid));
    CHECK_ERROR(cudaFree(beam_profile->radial_profile));
    CHECK_ERROR(cudaFree(beam_profile->cross_grid));
    CHECK_ERROR(cudaFree(beam_profile->cross_profile));
}

void free_global(TaigaGlobals *taiga_global){
    CHECK_ERROR(cudaFree(taiga_global->rad));
    CHECK_ERROR(cudaFree(taiga_global->z));
    CHECK_ERROR(cudaFree(taiga_global->tor));
    CHECK_ERROR(cudaFree(taiga_global->vrad));
    CHECK_ERROR(cudaFree(taiga_global->vz));
    CHECK_ERROR(cudaFree(taiga_global->vtor));
    CHECK_ERROR(cudaFree(taiga_global->detcellid));
    CHECK_ERROR(cudaFree(taiga_global->intensity));
    CHECK_ERROR(cudaFree(taiga_global->time_of_flight));
}

void free_common(TaigaCommons *shared_common){
    CHECK_ERROR(cudaFree(shared_common->grid_size));
    CHECK_ERROR(cudaFree(shared_common->spline_rgrid));
    CHECK_ERROR(cudaFree(shared_common->spline_zgrid));
    CHECK_ERROR(cudaFree(shared_common->brad));
    CHECK_ERROR(cudaFree(shared_common->bz));
    CHECK_ERROR(cudaFree(shared_common->btor));
    CHECK_ERROR(cudaFree(shared_common->erad));
    CHECK_ERROR(cudaFree(shared_common->ez));
    CHECK_ERROR(cudaFree(shared_common->etor));
    CHECK_ERROR(cudaFree(shared_common->psi_n));
    CHECK_ERROR(cudaFree(shared_common->detector_geometry));
    CHECK_ERROR(cudaFree(shared_common->ts_psi));
    CHECK_ERROR(cudaFree(shared_common->ts_density));
    CHECK_ERROR(cudaFree(shared_common->ts_temperature));
}

