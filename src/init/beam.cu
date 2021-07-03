#include <cuda.h>
#include "dataio/beam.h"

void init_coords(BeamProp *beam, ShotProp *shot, RunProp *run, TaigaGlobals *g_host, TaigaGlobals *g_shared) {
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
