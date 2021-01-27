void init_taiga(TaigaGlobals *g, TaigaCommons *s);
void init_host(TaigaGlobals *g, TaigaCommons *s);
void init_grid(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared, TaigaCommons *s_device);
void init_coords(TaigaGlobals *g_host, TaigaGlobals *g_shared, TaigaGlobals *g_device, BeamProp beam, ShotProp shot, RunProp run);
void init_beam_profile(BeamProfile *dev_prof, ShotProp shot);
void coord_memcopy(TaigaGlobals *g_host, TaigaGlobals *g, TaigaCommons *s, BeamProp beam, ShotProp shot, RunProp run);
void init_device_structs(TaigaGlobals *g_host, TaigaGlobals *g, TaigaCommons *s_host, TaigaCommons *s, BeamProp beam, ShotProp shot, RunProp run);
