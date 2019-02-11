char* clean_string (char* str_in);
void init_taiga_props(char* par_name, char* par_value, shot_prop *shot, beam_prop *beam);
int runnumber_reader(char* filename, shot_prop *shot);
int parameter_reader(char* filename, shot_prop *shot, beam_prop *beam);
