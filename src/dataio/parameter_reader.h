char* clean_string (char* str_in);
void init_taiga_props(char* par_name, char* par_value, BeamProp *beam, ShotProp *shot, RunProp *run);
int runnumber_reader(ShotProp *shot, RunProp *run);
int parameter_reader(BeamProp *beam, ShotProp *shot, RunProp *run);
void set_solver(RunProp *run, char* solver);
double get_mass(char *s);
