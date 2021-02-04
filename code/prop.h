#ifndef PROP_H
#define PROP_H

#define UNDEFINED_FLOAT -1e99
#define CALCULATION_NOT_FINISHED -1
#define CALCULATION_FINISHED 0
#define OUT_OF_RANGE -1

struct TaigaGlobals{
    double *rad, *z, *tor, *vrad, *vz, *vtor;
    int *detcellid;
    long particle_number;
};

struct TaigaLocals{
    double coords[6];
    int detcellid;
};

struct TaigaCommons{
    long max_step_number;        // N_step
    long step_counter;
    double eperm;
    double timestep;
    int *grid_size; 
    double *spline_rgrid;
    double *spline_zgrid;
    double **brad, **bz, **btor;
    bool is_electric_field_on;
    double **erad, **ez, **etor;
    double *detector_geometry;
};

struct BeamProp{
    char matter[40];
    double mass;                // in amu
    double energy;              // in keV
    double diameter;            // in meter
    double toroidal_deflection; // in radian
    double vertical_deflection; // in radian
    double deflection_coordinate;// radial position of deflection plates in meter
};

void init_beam_prop(BeamProp *beam){
    strcpy(beam->matter, "Li");
    beam->mass =  7.016004558;
    beam->toroidal_deflection = 0;
    beam->vertical_deflection = 0;
    beam->deflection_coordinate = 2.3;
}

struct ShotProp{
    char name[100];
    char shotnumber[40];
    char time[40];
    char detector_mask[40];
    char detector_geometry[100];
    bool is_electric_field_on;
};

void init_shot_prop(ShotProp *shot){
    strcpy(shot->name, "11774_1000");
    strcpy(shot->shotnumber, "11774");
    strcpy(shot->time, "1000");
    strcpy(shot->detector_mask, "test");
    strcpy(shot->detector_geometry, "0.685,0.23,0,38,0");
    shot->is_electric_field_on = false;
};

struct RunProp{
    int debug;
    int help;
    long particle_number;
    int block_number;
    int block_size;         //size of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)
    int step_host;          // on HDD
    int step_device;        // on GPU
    double timestep;
    double cpu_time_copy, cuda_time_copy, cuda_time_core;
    char runnumber[10];
    char parameter_file[200];
    char runnumber_file[200];
    char ion_source_file[200];
    char io_coordinate_order[200];
    char folder_out[200];
};

void init_run_prop(RunProp *run){
    run->debug = 0;
    run->help = 0;
    run->particle_number = 1;
    run->block_number = 1;
    run->block_size = 192;   //size of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)
    run->step_host = 1;      // on HDD
    run->step_device = 2000; // on GPU
    run->timestep = 1e-9;
    run->cpu_time_copy = UNDEFINED_FLOAT;
    run->cuda_time_copy = UNDEFINED_FLOAT;
    run->cuda_time_core = UNDEFINED_FLOAT;
    strcpy(run->runnumber, "0");
    strcpy(run->parameter_file, "parameters.sh");
    strcpy(run->runnumber_file, "runnumber");
    strcpy(run->ion_source_file, "");
    strcpy(run->io_coordinate_order, "rzt");
    strcpy(run->folder_out, "results/00000/");
};

struct DetectorProp{
    int detector_module_on;
    int length_xgrid;
    int length_ygrid;
    int number_of_detector_cells;
    int* counter;
    double* xgrid;
    double* ygrid;
};

struct BeamProfile{
    int radial_length;
    double *radial_grid;
    double *radial_profile;
    
    int cross_length;
    double *cross_grid;
    double *cross_profile;
};

#endif
