struct beam_prop{
    char matter[40];
    double mass;                // in amu
    double energy;              // in keV
    double diameter;            // in meter
    double toroidal_deflection; // in radian
    double vertical_deflection; // in radian    
};

void init_beam_prop(beam_prop beam){
    strcpy(beam.matter, "Li");
    beam.mass =  7.016004558;
    beam.energy = 60;
    beam.diameter = 0.025;
    beam.toroidal_deflection = 0;
    beam.vertical_deflection = 0;
}

struct shot_prop{
    char name[100];
    char shotnumber[40];
    char time[40];
    char detector_mask[40];
    char detector_geometry[100];
    int electric_field_module;
};

void init_shot_prop(shot_prop shot){
    strcpy(shot.name, "11774_1000");
    strcpy(shot.shotnumber, "11774");
    strcpy(shot.time, "1000");
    strcpy(shot.detector_mask, "test");
    strcpy(shot.detector_geometry, "0.685,0.23,0,38,0");
    shot.electric_field_module = 0;
};

struct run_prop{
    int runnumber;
    int debug;
    int help;
    int particle_number;
    int block_number;
    int block_size;         //size of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)
    int step_host;          // on HDD
    int step_device;        // on GPU
    double timestep;
    char parameter_file[200];
    char runnumber_file[200];
    char ion_source_file[200];
    char io_coordinate_order[200];
};

void init_run_prop(run_prop run){
    run.runnumber = 0;
    run.debug = 0;
    run.help = 0;
    run.particle_number = 1;
    run.block_number = 1;
    run.block_size = 192;   //size of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)
    run.step_host = 1;      // on HDD
    run.step_device = 2000; // on GPU
    run.timestep = 1e-9;
    strcpy(run.parameter_file, "parameters.sh");
    strcpy(run.runnumber_file, "runnumber");
    strcpy(run.ion_source_file, "");
    strcpy(run.io_coordinate_order, "rzt");
};


