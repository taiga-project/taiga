#include <cuda.h>
#include "taiga_init.h"
#include "taiga_constants.h"
#include "basic_functions.c"
#include "running/generate_coords.cuh"
#include "dataio/beam.h"

void init_taiga(TaigaGlobals *g, TaigaCommons *s){
}

void init_host(TaigaGlobals *g, TaigaCommons *s){
    s->step_counter = 0;
    s->max_step_number = 0;
    s->eperm = 0;
    s->timestep = 0;
    s->espline_on = false;
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
    
//    cudaMemcpy(s_device, s_shared, size_commons, cudaMemcpyHostToDevice);
    
    printf(" GRID SIZE: %d %d \n", s_host->grid_size[0], s_host->grid_size[1]);
}

void init_coords(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_host, TaigaGlobals *g_shared){
    size_t size_coord = run.block_size * run.block_number * sizeof(double);
    size_t size_detcellid = run.block_size * run.block_number * sizeof(int);
    size_t size_globals = sizeof(TaigaGlobals);
    
    if (!FASTMODE){
        double* shared_rad;
        double* shared_z;
        double* shared_tor;
        double* shared_vrad;
        double* shared_vz;
        double* shared_vtor;
        int* shared_detcellid;
        
        g_host->rad = (double*)malloc(size_coord);
        g_host->z   = (double*)malloc(size_coord);
        g_host->tor = (double*)malloc(size_coord);
        g_host->vrad = (double*)malloc(size_coord);
        g_host->vz   = (double*)malloc(size_coord);
        g_host->vtor = (double*)malloc(size_coord);
        g_host->detcellid = (int*)malloc(size_detcellid);
        load_beam(g_host, beam, shot, run);
        
        for (int i=0; i<run.block_size * run.block_number; ++i){
            g_host->detcellid[i] = -1;
        }
        
        memcpy(g_shared, g_host, size_globals);
        
        cudaMalloc((void **) &shared_rad, size_coord);
        cudaMalloc((void **) &shared_z,   size_coord);
        cudaMalloc((void **) &shared_tor, size_coord);
        cudaMalloc((void **) &shared_vrad, size_coord);
        cudaMalloc((void **) &shared_vz,   size_coord);
        cudaMalloc((void **) &shared_vtor, size_coord);
        cudaMalloc((void **) &shared_detcellid, size_detcellid);
        
        cudaMemcpy(shared_rad,       g_host->rad,       size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_z,         g_host->z,         size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_tor,       g_host->tor,       size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_vrad,      g_host->vrad,      size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_vz,        g_host->vz,        size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_vtor,      g_host->vtor,      size_coord,  cudaMemcpyHostToDevice);
        cudaMemcpy(shared_detcellid, g_host->detcellid, size_detcellid, cudaMemcpyHostToDevice);
        
        g_shared->rad  = shared_rad;
        g_shared->z    = shared_z;
        g_shared->tor  = shared_tor;
        g_shared->vrad = shared_vrad;
        g_shared->vz   = shared_vz;
        g_shared->vtor = shared_vtor;
        g_shared->detcellid = shared_detcellid;
    }else{
        BeamProfile dev_beam_prof;
        init_beam_profile(&dev_beam_prof, shot);
        printf("i84 \n");//# 
        //generate_coords <<< run.block_number, run.block_size >>> (*g, *s, beam, dev_beam_prof);
    }
}

void init_beam_profile(BeamProfile *dev_prof, ShotProp shot){
    /*#BeamProfile host_prof;
    
//#    init_ion_profile(shot.name, host_prof);

    size_t dimBF = sizeof(beam_distribution);
    size_t dimBR = sizeof(double)*host_prof.radial.N;
    size_t dimBX = sizeof(double)*host_prof.cross_section.N;

//#    load_ion_profile(shot.name, host_prof);
printf("i87 \n");
    cudaMalloc((void **) &(dev_prof->radial), dimBF);
    cudaMalloc((void **) &(dev_prof->radial.grid),    dimBR);
    cudaMalloc((void **) &(dev_prof->radial.profile), dimBR);
    cudaMemcpy(&(dev_prof->radial.grid),    &(host_prof.radial.grid),    dimBR, cudaMemcpyHostToDevice);
    cudaMemcpy(&(dev_prof->radial.profile), &(host_prof.radial.profile), dimBR, cudaMemcpyHostToDevice);
printf("i93 \n");
    cudaMalloc((void **) &(dev_prof->cross_section), dimBF);
    cudaMalloc((void **) &(dev_prof->cross_section.grid),    dimBX);
    cudaMalloc((void **) &(dev_prof->cross_section.profile), dimBX);
    cudaMemcpy(&(dev_prof->cross_section.grid),    &(host_prof.cross_section.grid),    dimBX, cudaMemcpyHostToDevice);
    cudaMemcpy(&(dev_prof->cross_section.profile), &(host_prof.cross_section.profile), dimBX, cudaMemcpyHostToDevice);*/ 
}

void coord_memcopy_back(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_shared, TaigaGlobals *g_host, TaigaCommons *c_shared, TaigaCommons *c_host){
    size_t size_coord = run.block_size * run.block_number * sizeof(double);
    size_t size_detcellid = run.block_size * run.block_number * sizeof(int);
    cudaMemcpy(g_shared->rad, g_host->rad, size_coord, cudaMemcpyDeviceToHost);
    cudaMemcpy(g_shared->z,   g_host->z,   size_coord, cudaMemcpyDeviceToHost);
    cudaMemcpy(g_shared->tor, g_host->tor, size_coord, cudaMemcpyDeviceToHost);
    cudaMemcpy(g_shared->vrad, g_host->vrad, size_coord, cudaMemcpyDeviceToHost);
    cudaMemcpy(g_shared->vz,   g_host->vz,   size_coord, cudaMemcpyDeviceToHost);
    cudaMemcpy(g_shared->vtor, g_host->vtor, size_coord, cudaMemcpyDeviceToHost);
    cudaMemcpy(g_shared->detcellid, g_host->detcellid, size_detcellid, cudaMemcpyDeviceToHost);
    /*cudaMemcpy(&(g_host->rad),       &(g->rad),       size_coord,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->z),         &(g->z),         size_coord,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->tor),       &(g->tor),       size_coord,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->vrad),      &(g->vrad),      size_coord,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->vz),        &(g->vz),        size_coord,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->vtor),      &(g->vtor),      size_coord,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->detcellid), &(g->detcellid), size_detcellid, cudaMemcpyHostToDevice);*/
    
}

void init_device_structs(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *g_device, TaigaGlobals *g_shared, TaigaCommons *c_device, TaigaCommons *c_shared){
printf("i143");
    g_shared->particle_number = run.particle_number;
    c_shared->max_step_number = run.step_device;        // N_step
    c_shared->step_counter    = 0;
    c_shared->eperm           = ELEMENTARY_CHARGE/ AMU/ beam.mass;
    c_shared->timestep        = run.timestep;
//    cudaMemcpy(g_device, g_shared, sizeof(TaigaGlobals), cudaMemcpyHostToDevice);
//    cudaMemcpy(c_device, c_shared, sizeof(TaigaCommons), cudaMemcpyHostToDevice);
printf("i155");
}

void set_particle_number(TaigaGlobals *host_global, RunProp *run){
    if (READINPUTPROF == 1){
        double *X_temp;
        host_global->particle_number = read_vector(&X_temp, "input", "manual_profile", "rad.dat");
        run->block_number = host_global->particle_number / run->block_size+1;
        free(X_temp);
    }else{
        host_global->particle_number = run->block_size * run->block_number;
    }
}
