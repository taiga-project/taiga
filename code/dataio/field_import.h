int spline_read_and_init(ShotProp shot, RunProp run, char* field_name, double ***return_s_ptr, int dimRZ);
int magnetic_field_read_and_init(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared);
int electric_field_read_and_init(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared);
