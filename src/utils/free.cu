#include "utils/free.cuh"

void free_taiga() {

}


void free_detector(DetectorProp *detector, int where){
    FREE(where, detector->counter);
    FREE(where, detector->xgrid);
    FREE(where, detector->ygrid);
}

void free_beam(BeamProfile *beam_profile, int where){
    FREE(where, beam_profile->radial_grid);
    FREE(where, beam_profile->radial_profile);
    FREE(where, beam_profile->cross_grid);
    FREE(where, beam_profile->cross_profile);
}

void free_global(TaigaGlobals *taiga_global, int where){
    FREE(where, taiga_global->rad);
    FREE(where, taiga_global->z);
    FREE(where, taiga_global->tor);
    FREE(where, taiga_global->vrad);
    FREE(where, taiga_global->vz);
    FREE(where, taiga_global->vtor);
    FREE(where, taiga_global->detcellid);
    FREE(where, taiga_global->intensity);
    FREE(where, taiga_global->time_of_flight);
}

void free_common(TaigaCommons *shared_common, int where){
    FREE(where, shared_common->grid_size);
    FREE(where, shared_common->spline_rgrid);
    FREE(where, shared_common->spline_zgrid);
    FREE(where, shared_common->brad);
    FREE(where, shared_common->bz);
    FREE(where, shared_common->btor);
    FREE(where, shared_common->erad);
    FREE(where, shared_common->ez);
    FREE(where, shared_common->etor);
    FREE(where, shared_common->psi_n);
    FREE(where, shared_common->detector_geometry);
    FREE(where, shared_common->ts_psi);
    FREE(where, shared_common->ts_density);
    FREE(where, shared_common->ts_temperature);
}

