#ifndef PROP_H
#define PROP_H

typedef struct TaigaGlobals{
    double *rad, *z, *tor, *vrad, *vz, *vtor;
    int *detcellid;
    long particle_number;
}TaigaGlobals;

typedef struct TaigaCommons{
    long max_step_number;        // N_step
    long step_counter;
    double eperm;
    double timestep;
    int solver;
    int field_interpolation_method;
    int *grid_size;
    double *spline_rgrid;
    double *spline_zgrid;
    double **brad, **bz, **btor;
    bool is_electric_field_on;
    double **erad, **ez, **etor;
    bool is_magnetic_field_perturbation;
    double **psi_n;
    double *detector_geometry;
}TaigaCommons;

typedef struct BeamProp{
    char matter[STRING_LENGTH];
    double mass;                // in amu
    double energy;              // in keV
    double diameter;            // in meter
    double toroidal_deflection; // in radian
    double vertical_deflection; // in radian
    double deflection_coordinate;// radial position of deflection plates in meter
}BeamProp;

typedef struct ShotProp{
    char name[STRING_LENGTH];
    char shotnumber[STRING_LENGTH];
    char time[STRING_LENGTH];
    char detector_mask[STRING_LENGTH];
    char detector_geometry[STRING_LENGTH];
}ShotProp;

typedef struct RunProp{
    int debug;
    int help;
    long particle_number;
    int block_number;
    int block_size;         //size of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)
    long step_host;          // on HDD
    long step_device;        // on GPU
    double timestep;
    int solver;
    int field_interpolation_method;
    bool is_electric_field_on;
    bool is_magnetic_field_perturbation;
    double cpu_time_copy, cuda_time_copy, cuda_time_core;
    char runnumber[STRING_LENGTH];
    char parameter_file[STRING_LENGTH];
    char runnumber_file[STRING_LENGTH];
    char ion_source_file[STRING_LENGTH];
    char io_coordinate_order[STRING_LENGTH];
    char folder_out[STRING_LENGTH];
}RunProp;

typedef struct DetectorProp{
    int detector_module_on;
    int length_xgrid;
    int length_ygrid;
    int number_of_detector_cells;
    long* counter;
    double* xgrid;
    double* ygrid;
}DetectorProp;

typedef struct BeamProfile{
    long radial_length;
    double *radial_grid;
    double *radial_profile;
    long cross_length;
    double *cross_grid;
    double *cross_profile;
}BeamProfile;

#endif
