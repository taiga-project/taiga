#include "utils/prop.h"

#include <string.h>

void init_beam_prop(BeamProp *beam){
    strcpy(beam->species, "Li");
    beam->charge = 1.0;
    beam->toroidal_deflection = 0.0;
    beam->vertical_deflection = 0.0;
    beam->deflection_radial_coordinate = 2.3;
}

void init_shot_prop(ShotProp *shot){
    strcpy(shot->name, "11774_1000");
    strcpy(shot->shotnumber, "11774");
    strcpy(shot->time, "1000");
    strcpy(shot->detector_mask, "test");
    strcpy(shot->detector_geometry, "0.685,0.23,0,38,0");
}

void init_run_prop(RunProp *run){
    run->mode = ALL_IO;
    run->debug = 0;
    run->help = 0;
    run->particle_number = 1;
    run->block_number = 1;
    run->block_size = 192;   //size of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)
    run->step_host = 1;      // on HDD
    run->step_device = 2000; // on GPU
    run->timestep = 1e-9;
    run->solver = SOLVER_RK45;
    run->init_source = READ_COORDINATES;
    run->field_interpolation_method = CUBIC_SPLINE;
    run->detect_interpolation_method = LINEAR_INTERPOLATION;
    run->is_electric_field_on = false;
    run->is_magnetic_field_perturbation = false;
    run->is_ionisation_on = false;
    run->cpu_time_copy = UNDEFINED_FLOAT;
    run->cuda_time_copy = UNDEFINED_FLOAT;
    run->cuda_time_core = UNDEFINED_FLOAT;
    strcpy(run->runnumber, UNDEFINED_RUNNUMBER);
    strcpy(run->parameter_file, "parameters.sh");
    strcpy(run->runnumber_file, "runnumber");
    strcpy(run->ion_source_file, "");
    strcpy(run->io_coordinate_order, "rzt");
    strcpy(run->folder_out, "results/00000/");
}