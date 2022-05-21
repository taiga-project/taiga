#ifndef TAIGA_FREE_CUH
#define TAIGA_FREE_CUH

#define HOST 0
#define SHARED 1
#define DEVICE 2
#define FREE(where, p) { \
    if (where == HOST) {free(p);} \
    else {CHECK_ERROR(cudaFree(p));} \
}
#define FREE_MAIN(where, p) { \
    if (where == DEVICE) {CHECK_ERROR(cudaFree(p));} \
    else {free(p);} \
}

void free_taiga(TaigaGlobals *host_global, TaigaGlobals *shared_global, TaigaGlobals *device_global,
                TaigaCommons *host_common,  TaigaCommons *shared_common,  TaigaCommons *device_common,
                DetectorProp *shared_detector, DetectorProp *device_detector,
                ShotProp *shot, BeamProp *beam, RunProp *run,
                double *host_service_array, double *device_service_array);

void free_detector(DetectorProp *detector, int where);
void free_beam(BeamProfile *beam_profile, int where);
void free_global(TaigaGlobals *taiga_global, int where);
void free_common(TaigaCommons *shared_common, int where);

#endif //TAIGA_FREE_CUH
