#include <cuda.h>
#include "taiga_init.h"
#include "taiga_constants.h"
#include "basic_functions.c"
#include "running/generate_coords.cuh"
#include "dataio/beam.h"

void init_host(TaigaGlobals *g, TaigaCommons *s){
    s->step_counter = 0;
    s->max_step_number = 0;
    s->eperm = 0;
    s->timestep = 0;
    s->is_electric_field_on = false;
    g->particle_number = 0;
}

void init_grid(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared){
    
    int* shared_grid_size;
    double *shared_rgrid, *shared_zgrid;
    size_t size_commons = sizeof(TaigaCommons);
    size_t size_grid_dim = 2*sizeof(int);
    
    s_host->grid_size = (int*)malloc(size_grid_dim);
    s_host->grid_size[0] = read_vector(&s_host->spline_rgrid, "input/fieldSpl", shot.name, "r.spline");
    size_t size_R = s_host->grid_size[0] * sizeof(double);
    s_host->grid_size[1] = read_vector(&s_host->spline_zgrid, "input/fieldSpl", shot.name, "z.spline");
    size_t size_Z = s_host->grid_size[1] * sizeof(double);
    
    memcpy(s_shared, s_host, size_commons);
    
    cudaMalloc((void **) &shared_grid_size, size_grid_dim);
    cudaMemcpy(shared_grid_size, s_host->grid_size, size_grid_dim, cudaMemcpyHostToDevice);
    s_shared->grid_size = shared_grid_size;
    
    cudaMalloc((void **) &shared_rgrid, size_R);
    cudaMemcpy(shared_rgrid, s_host->spline_rgrid, size_R, cudaMemcpyHostToDevice);
    s_shared->spline_rgrid = shared_rgrid;
    
    cudaMalloc((void **) &shared_zgrid, size_Z);
    cudaMemcpy(shared_zgrid, s_host->spline_zgrid, size_Z, cudaMemcpyHostToDevice);
    s_shared->spline_zgrid = shared_zgrid;
    
    printf(" GRID SIZE: %d %d \n", s_host->grid_size[0], s_host->grid_size[1]);
}

void init_coords(BeamProp *beam, ShotProp *shot, RunProp *run, TaigaGlobals *g_host, TaigaGlobals *g_shared){
    size_t size_coord = run->block_size * run->block_number * sizeof(double);
    size_t size_detcellid = run->block_size * run->block_number * sizeof(int);
    size_t size_globals = sizeof(TaigaGlobals);
    
    double* shared_rad;
    double* shared_z;
    double* shared_tor;
    double* shared_vrad;
    double* shared_vz;
    double* shared_vtor;
    int* shared_detcellid;
    
    if (!FASTMODE){
        g_host->rad = (double*)malloc(size_coord);
        g_host->z   = (double*)malloc(size_coord);
        g_host->tor = (double*)malloc(size_coord);
        g_host->vrad = (double*)malloc(size_coord);
        g_host->vz   = (double*)malloc(size_coord);
        g_host->vtor = (double*)malloc(size_coord);
        g_host->detcellid = (int*)malloc(size_detcellid);
        load_beam(g_host, beam, shot, run);
        
        for (int i=0; i<run->block_size * run->block_number; ++i){
            g_host->detcellid[i] = CALCULATION_NOT_FINISHED;
        }
        
        memcpy(g_shared, g_host, size_globals);
    }
    
    cudaMalloc((void **) &shared_rad, size_coord);
    cudaMalloc((void **) &shared_z,   size_coord);
    cudaMalloc((void **) &shared_tor, size_coord);
    cudaMalloc((void **) &shared_vrad, size_coord);
    cudaMalloc((void **) &shared_vz,   size_coord);
    cudaMalloc((void **) &shared_vtor, size_coord);
    cudaMalloc((void **) &shared_detcellid, size_detcellid);
    
    if (!FASTMODE){
        cudaMemcpy(shared_rad,       g_host->rad,       size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_z,         g_host->z,         size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_tor,       g_host->tor,       size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_vrad,      g_host->vrad,      size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_vz,        g_host->vz,        size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_vtor,      g_host->vtor,      size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_detcellid, g_host->detcellid, size_detcellid, cudaMemcpyHostToDevice);
    }
    
    g_shared->rad  = shared_rad;
    g_shared->z    = shared_z;
    g_shared->tor  = shared_tor;
    g_shared->vrad = shared_vrad;
    g_shared->vz   = shared_vz;
    g_shared->vtor = shared_vtor;
    g_shared->detcellid = shared_detcellid;
}

void init_fastmode(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *device_global){
    BeamProfile *device_prof;
    size_t size_prof = sizeof(BeamProfile);
    cudaMalloc((void **) &device_prof, size_prof);
    init_beam_profile(device_prof, shot);
    generate_coords <<< run.block_number, run.block_size >>> (device_global, beam, device_prof);
}

void init_beam_profile(BeamProfile *device_prof, ShotProp shot){
    BeamProfile *host_prof, *shared_prof;
    size_t size_prof = sizeof(BeamProfile);
    host_prof = (BeamProfile*)malloc(size_prof);
    shared_prof = (BeamProfile*)malloc(size_prof);
    init_ion_profile(shot.name, host_prof);
    size_t size_rad_prof = sizeof(double)*host_prof->radial_length;
    size_t size_cross_prof = sizeof(double)*host_prof->cross_length;
    
    double *shared_radial_grid;
    double *shared_radial_profile;
    double *shared_cross_grid;
    double *shared_cross_profile;
    
    cudaMalloc((void **) &shared_radial_grid,    size_rad_prof);
    cudaMalloc((void **) &shared_radial_profile, size_rad_prof);
    cudaMalloc((void **) &shared_cross_grid,     size_cross_prof);
    cudaMalloc((void **) &shared_cross_profile,  size_cross_prof);
    
    cudaMemcpy(shared_radial_grid,    host_prof->radial_grid,    size_rad_prof, cudaMemcpyHostToDevice);
    cudaMemcpy(shared_radial_profile, host_prof->radial_profile, size_rad_prof, cudaMemcpyHostToDevice);
    cudaMemcpy(shared_cross_grid,     host_prof->cross_grid,     size_cross_prof, cudaMemcpyHostToDevice);
    cudaMemcpy(shared_cross_profile,  host_prof->cross_profile,  size_cross_prof, cudaMemcpyHostToDevice);
    
    shared_prof->radial_grid    = shared_radial_grid;
    shared_prof->radial_profile = shared_radial_profile;
    shared_prof->cross_grid     = shared_cross_grid;
    shared_prof->cross_profile  = shared_cross_profile;
    
    cudaMemcpy(device_prof,  shared_prof,  size_prof, cudaMemcpyHostToDevice);
}

void coord_memcopy_back(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_host, TaigaGlobals *g_shared){
    size_t size_coord = run.block_size * run.block_number * sizeof(double);
    size_t size_detcellid = run.block_size * run.block_number * sizeof(int);
    
    double* host_rad =(double*)malloc(size_coord);
    double* host_z =(double*)malloc(size_coord);
    double* host_tor =(double*)malloc(size_coord);
    double* host_vrad =(double*)malloc(size_coord);
    double* host_vz =(double*)malloc(size_coord);
    double* host_vtor =(double*)malloc(size_coord);
    int* host_detcellid =(int*)malloc(size_detcellid);
    
    cudaMemcpy(host_rad,  g_shared->rad,  size_coord, cudaMemcpyDeviceToHost); g_host->rad = host_rad;
    cudaMemcpy(host_z,    g_shared->z,    size_coord, cudaMemcpyDeviceToHost); g_host->z = host_z;
    cudaMemcpy(host_tor,  g_shared->tor,  size_coord, cudaMemcpyDeviceToHost); g_host->tor = host_tor;
    cudaMemcpy(host_vrad, g_shared->vrad, size_coord, cudaMemcpyDeviceToHost); g_host->vrad = host_vrad;
    cudaMemcpy(host_vz,   g_shared->vz,   size_coord, cudaMemcpyDeviceToHost); g_host->vz = host_vz;
    cudaMemcpy(host_vtor, g_shared->vtor, size_coord, cudaMemcpyDeviceToHost); g_host->vtor = host_vtor;
    cudaMemcpy(host_detcellid, g_shared->detcellid, size_detcellid, cudaMemcpyDeviceToHost); g_host->detcellid = host_detcellid;
}

void sync_device_structs(TaigaGlobals *g_device, TaigaGlobals *g_shared, TaigaCommons *c_device, TaigaCommons *c_shared){
    cudaMemcpy(c_device, c_shared, sizeof(TaigaCommons), cudaMemcpyHostToDevice);
    if (!FASTMODE){
        cudaMemcpy(g_device, g_shared, sizeof(TaigaGlobals), cudaMemcpyHostToDevice);
    }
}

void init_device_structs(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_shared, TaigaCommons *c_shared){
    g_shared->particle_number       = run.particle_number;
    c_shared->max_step_number       = run.step_device;
    c_shared->step_counter          = 0;
    c_shared->eperm                 = ELEMENTARY_CHARGE/ AMU/ beam.mass;
    c_shared->timestep              = run.timestep;
    c_shared->solver                = run.solver;
}

void set_particle_number(RunProp *run, TaigaGlobals *host_global, TaigaGlobals *shared_global){
    if (READINPUTPROF == 1){
        double *X_temp;
        host_global->particle_number = read_vector(&X_temp, "input", "manual_profile", "rad.dat");
        run->block_number = host_global->particle_number / run->block_size+1;
        free(X_temp);
    }else{
        host_global->particle_number = run->block_size * run->block_number;
    }
    shared_global->particle_number = host_global->particle_number;
}

void set_detector_geometry(ShotProp shot, TaigaCommons *host_common, TaigaCommons *shared_common){
    char *tokaniser;
    double *shared_detector_geometry;
    
    size_t size_detector = 5 * sizeof(double);
    host_common->detector_geometry = (double *)malloc(size_detector);
    tokaniser = strtok(shot.detector_geometry, ",");
        host_common->detector_geometry[0] = strtod (tokaniser, NULL);
    tokaniser = strtok(NULL,",");
        host_common->detector_geometry[1] = strtod (tokaniser, NULL);
    tokaniser = strtok(NULL,",");
        host_common->detector_geometry[2] = strtod (tokaniser, NULL);
    tokaniser = strtok(NULL,",");
        host_common->detector_geometry[3] = (strtod (tokaniser, NULL) * PI/180.0);
    tokaniser = strtok(NULL,",");
        host_common->detector_geometry[4] = (strtod (tokaniser, NULL) * PI/180.0);
    
    cudaMalloc((void **) &shared_detector_geometry, size_detector);
    cudaMemcpy(shared_detector_geometry, host_common->detector_geometry, size_detector, cudaMemcpyHostToDevice);
    shared_common->detector_geometry = shared_detector_geometry;
}

void init_detector(DetectorProp* shared_detector, DetectorProp *device_detector, ShotProp shot){
    double *host_detector_xgrid, *device_detector_xgrid;
    double *host_detector_ygrid, *device_detector_ygrid;
    long *device_counter;
    
    shared_detector->length_xgrid = read_vector(&host_detector_xgrid, "input/detector", shot.detector_mask, "detx", false)/2;
    shared_detector->length_ygrid = read_vector(&host_detector_ygrid, "input/detector", shot.detector_mask, "dety", false)/2;
    shared_detector->number_of_detector_cells = shared_detector->length_xgrid * shared_detector->length_ygrid;

    shared_detector->detector_module_on = (( shared_detector->length_xgrid > 0) && ( shared_detector->length_ygrid > 0));
    if (shared_detector->detector_module_on) {
        printf("===============================\n");
        printf("Detector postprocessor module: ON (%d x %d / %s)\n", shared_detector->length_xgrid, shared_detector->length_ygrid, shot.detector_mask);
        size_t size_detector_xgrid = 2 * shared_detector->length_xgrid * sizeof(double);
        size_t size_detector_ygrid = 2 * shared_detector->length_ygrid * sizeof(double);
        size_t size_counter = shared_detector->number_of_detector_cells * sizeof(long);
        size_t size_detector = sizeof(DetectorProp);
        
        cudaMalloc((void **) &device_detector_xgrid, size_detector_xgrid);
        cudaMalloc((void **) &device_detector_ygrid, size_detector_ygrid);
        cudaMalloc((void **) &device_counter, size_counter);
        
        cudaMemcpy(device_detector_xgrid, host_detector_xgrid, size_detector_xgrid, cudaMemcpyHostToDevice);
        cudaMemcpy(device_detector_ygrid, host_detector_ygrid, size_detector_ygrid, cudaMemcpyHostToDevice);
        shared_detector->xgrid = device_detector_xgrid;
        shared_detector->ygrid = device_detector_ygrid;
        shared_detector->counter = device_counter;
        
        cudaMemcpy(device_detector, shared_detector, size_detector, cudaMemcpyHostToDevice);
    }else{
        printf("===============================\n");
        printf("Detector postprocessor module: OFF\n");
    }
}
