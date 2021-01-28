void init_taiga(TaigaGlobals *g, TaigaCommons *s);
void init_host(TaigaGlobals *g, TaigaCommons *s);
void init_grid(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared);
void init_coords(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_host, TaigaGlobals *g_shared);
void init_beam_profile(BeamProfile *dev_prof, ShotProp shot);
void coord_memcopy_back(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_shared, TaigaGlobals *g_host, TaigaCommons *c_shared, TaigaCommons *c_host);
void init_device_structs(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_device, TaigaGlobals *g_shared, TaigaCommons *c_device, TaigaCommons *c_shared);
void set_particle_number(TaigaGlobals *host_global, RunProp *run);
