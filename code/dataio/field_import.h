int spline_read_and_init(shot_prop shot, run_prop run, char* field_name, double ***return_s_ptr, int dimRZ);
int magnetic_field_read_and_init(shot_prop shot, run_prop run, double ***return_br_ptr, double ***return_bz_ptr, double ***return_bt_ptr, int dimRZ);
int electric_field_read_and_init(shot_prop shot, run_prop run, double ***return_er_ptr, double ***return_ez_ptr, double ***return_et_ptr, int dimRZ);
