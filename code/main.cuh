struct beam_prop{
    char matter[40] = "Li";
    double mass = 7.016004558;      // in amu
    double energy = 60;             // in keV
    double diameter = 0.025;        // in meter
    double toroidal_deflection = 0; // in radian
    double vertical_deflection = 0; // in radian    
};

struct shot_prop{
    char name[100] = "11774_1000";
    char shotnumber[40] = "11774";
    char time[40] = "1000";
    char detector_mask[40] = "test";
    char detector_geometry[100] = "0.685,0.23,0,38,0";
    int electric_field_module = 0;
};

struct run_prop{
    int runnumber = 0;
    int debug = 0;
    int help = 0;
    int particle_number = 1;
    int block_number = 1;
    int block_size = 192;   //size of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)
    int step_host = 1;      // on HDD
    int step_device = 2000; // on GPU
    char parameter_file[200] = "parameters.sh";
    char runnumber_file[200] = "runnumber";
};

void set_cuda(int debug_flag);
inline void cErrorCheck(const char *file, int line);

int main(int argc, char *argv[]);
void input_init_taiga(int argc, char *argv[], shot_prop *shot, beam_prop *beam, run_prop *run);

double get_mass(char *s);
void set_detector_geometry(double *DETECTOR, char* values);
void detector_module(double **x_ptr, double *detector, int *detcellid, char *detector_name, int max_blocks, int shot_block_size, int number_of_particles, char *export_folder, char *runnumber);

int spline_read_and_init(shot_prop shot, run_prop run, char* field_name, double ***return_s_ptr, int dimRZ);

int magnetic_field_read_and_init(shot_prop shot, run_prop run, double ***return_br_ptr, double ***return_bz_ptr, double ***return_bt_ptr, int dimRZ);
int electric_field_read_and_init(shot_prop shot, run_prop run, double ***return_er_ptr, double ***return_ez_ptr, double ***return_et_ptr, int dimRZ);

void print_help_message();
