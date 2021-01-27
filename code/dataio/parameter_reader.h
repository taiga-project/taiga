char* clean_string (char* str_in);
void init_taiga_props(char* par_name, char* par_value, ShotProp *shot, BeamProp *beam, RunProp *run);
int runnumber_reader(ShotProp *shot, RunProp *run);
int parameter_reader(ShotProp *shot, BeamProp *beam, RunProp *run);

double get_mass(char *s);
void set_detector_geometry(double *DETECTOR, char* values);
