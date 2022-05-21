#ifndef TAIGA_FREE_CUH
#define TAIGA_FREE_CUH

#define HOST 0
#define DEVICE 1
#define FREE(where, p) { \
    if (where == DEVICE) {CHECK_ERROR(cudaFree(p));} \
    else {free(p);} \
}

void free_detector(DetectorProp *detector, int where);
void free_beam(BeamProfile *beam_profile, int where);
void free_global(TaigaGlobals *taiga_global, int where);
void free_common(TaigaCommons *shared_common, int where);

#endif //TAIGA_FREE_CUH
