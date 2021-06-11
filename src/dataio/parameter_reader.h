char* clean_string (char* str);
void init_taiga_props(char* par_name, char* par_value, BeamProp *beam, ShotProp *shot, RunProp *run);
int runnumber_reader(ShotProp *shot, RunProp *run);
int parameter_reader(BeamProp *beam, ShotProp *shot, RunProp *run);
void set_solver(RunProp *run, char* solver);
void set_field_interpolation_method(RunProp *run, char* method);
double get_mass(char *name_of_ion);
