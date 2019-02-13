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

int main(int argc, char *argv[]);
void print_help_message();
void input_init_taiga(int argc, char *argv[], shot_prop *shot, beam_prop *beam, run_prop *run);
