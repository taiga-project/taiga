void init_taiga(device_global *g, device_shared *s);
void init_device(device_global *g, device_shared *s);
void init_grid(shot_prop shot, run_prop run, device_shared *s_host, device_shared *s);
void init_coords(device_global *g_host, device_global *g, device_shared *s, beam_prop beam, shot_prop shot, run_prop run);
void init_beam_profile(beam_profile *dev_prof, shot_prop shot);
void coord_memcopy(device_global *g_host, device_global *g, device_shared *s, beam_prop beam, shot_prop shot, run_prop run);
void init_device_structs(device_global *g_host, device_global *g, device_shared *s_host,  device_shared *s, beam_prop beam, shot_prop shot, run_prop run);
