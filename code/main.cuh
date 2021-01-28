int main(int argc, char *argv[]);
void print_help_message();
void fill_header_file(ShotProp shot, BeamProp beam, RunProp run, TaigaCommons *common);
void input_init_taiga(int argc, char *argv[], ShotProp *shot, BeamProp *beam, RunProp *run);
void save_trajectories(TaigaGlobals *host_global, RunProp run);
void save_endpoints(TaigaGlobals *host_global, RunProp run);
