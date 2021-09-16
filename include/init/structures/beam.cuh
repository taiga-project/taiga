#ifndef INIT_BEAM_CUH
#define INIT_BEAM_CUH

void init_coords(BeamProp *beam, ShotProp *shot, RunProp *run, TaigaGlobals *g_host, TaigaGlobals *g_shared);
void init_beam_profile(BeamProfile *device_prof, ShotProp shot);

#endif //INIT_BEAM_CUH