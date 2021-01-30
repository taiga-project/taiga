int main(int argc, char *argv[]);
void print_help_message();
void print_run_details(TaigaGlobals *host_global, TaigaCommons *host_common, ShotProp shot, RunProp run);
void fill_header_file(TaigaCommons *common, BeamProp beam, ShotProp shot, RunProp run);
void input_init_taiga(int argc, char *argv[], ShotProp *shot, BeamProp *beam, RunProp *run);
void save_trajectories(TaigaGlobals *host_global, RunProp run);
void save_endpoints(TaigaGlobals *host_global, RunProp run);
