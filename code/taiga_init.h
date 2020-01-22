void init_taiga(taiga_globals *g, taiga_commons *s);
void init_device(taiga_globals *g, taiga_commons *s);
void init_grid(shot_prop shot, run_prop run, taiga_commons *s_host, taiga_commons *s);
void init_coords(taiga_globals *g_host, taiga_globals *g, taiga_commons *s, beam_prop beam, shot_prop shot, run_prop run);
void init_beam_profile(beam_profile *dev_prof, shot_prop shot);
void coord_memcopy(taiga_globals *g_host, taiga_globals *g, taiga_commons *s, beam_prop beam, shot_prop shot, run_prop run);
void init_device_structs(taiga_globals *g_host, taiga_globals *g, taiga_commons *s_host,  taiga_commons *s, beam_prop beam, shot_prop shot, run_prop run);
