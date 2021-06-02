#ifndef INIT_H
#define INIT_H

void init_host(TaigaGlobals *g, TaigaCommons *s);
void init_grid(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared);
void init_device_structs(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_shared, TaigaCommons *c_shared);
void set_particle_number(RunProp *run, TaigaGlobals *host_global, TaigaGlobals *shared_global);

#endif //INIT_H
