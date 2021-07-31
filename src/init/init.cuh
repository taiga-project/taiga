#ifndef INIT_H
#define INIT_H

void init_host(TaigaGlobals *host_global, TaigaCommons *host_common);
void init_grid(ShotProp shot, RunProp run, TaigaCommons *host_common, TaigaCommons *shared_common);
void init_device_structs(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *shared_global, TaigaCommons *shared_common);
void set_particle_number(RunProp *run, TaigaGlobals *host_global, TaigaGlobals *shared_global);

#endif //INIT_H
