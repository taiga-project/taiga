#ifndef INIT_SYNC_CUH
#define INIT_SYNC_CUH

void sync_device_structs(TaigaGlobals *g_device, TaigaGlobals *g_shared, TaigaCommons *c_device, TaigaCommons *c_shared,
                         bool is_all_io);
void coord_memcopy_back(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_host, TaigaGlobals *g_shared);

#endif //INIT_SYNC_CUH
