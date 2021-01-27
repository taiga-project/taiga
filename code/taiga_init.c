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

void init_grid(ShotProp shot, RunProp run, TaigaCommons *s_host, TaigaCommons *s_shared, TaigaCommons *s_device){
    
    int* shared_grid_size;
    double *shared_rgrid, *shared_zgrid;
    size_t dimCommons = sizeof(TaigaCommons);
    size_t dimGP = 2*sizeof(int);
    
    
    printf("Reading init0\n");
    s_host->grid_size = (int*)malloc(dimGP);
    s_host->grid_size[0] = read_vector(&s_host->spline_rgrid, "input/fieldSpl", shot.name, "r.spline");
    size_t dimR = s_host->grid_size[0] * sizeof(double);
    s_host->grid_size[1] = read_vector(&s_host->spline_zgrid, "input/fieldSpl", shot.name, "z.spline");
    size_t dimZ = s_host->grid_size[1] * sizeof(double);
    
    printf("Reading init1\n");
    
    
    /*
    int *s_global__detcellid;
    cudaMalloc((void **) &(s_global__detcellid), dim_detcellid);
    
    h_global->particle_number = 10;
    h_global->detcellid[DETCELLID_INDEX] = 42;
    
    memcpy(s_global, h_global, dim_global);
    cudaMemcpy(s_global__detcellid, h_global->detcellid, dim_detcellid, cudaMemcpyHostToDevice);
    s_global->detcellid = s_global__detcellid;
    cudaMemcpy(d_global, s_global, dim_global, cudaMemcpyHostToDevice);
    */
    
    memcpy(s_shared, s_host, dimCommons);
    
    cudaMalloc((void **) &shared_grid_size, dimGP);
    cudaMemcpy(shared_grid_size, s_host->grid_size, dimGP, cudaMemcpyHostToDevice);
    s_shared->grid_size = shared_grid_size;
    
    cudaMalloc((void **) &shared_rgrid, dimR);
    cudaMemcpy(shared_rgrid, s_host->spline_rgrid, dimR, cudaMemcpyHostToDevice);
    s_shared->spline_rgrid = shared_rgrid;

    cudaMalloc((void **) &shared_zgrid, dimZ);
    cudaMemcpy(shared_zgrid, s_host->spline_zgrid, dimZ, cudaMemcpyHostToDevice);
    s_shared->spline_zgrid = shared_zgrid;
    
    //cudaMemcpy(&(s_device_grid_size),    &(s_host->grid_size),    dimGP, cudaMemcpyHostToDevice);
    //s_shared->

    cudaMemcpy(s_device, s_shared, dimCommons, cudaMemcpyHostToDevice);
    
    printf(" GRID SIZE: %d %d \n", s_host->grid_size[0], s_host->grid_size[1]);
}

void init_coords(TaigaGlobals *g_host, TaigaGlobals *g_shared, TaigaGlobals *g_device, BeamProp beam, ShotProp shot, RunProp run){
    size_t dimX = run.block_size * run.block_number * sizeof(double);
    size_t dimXI = run.block_size * run.block_number * sizeof(int);
    double* g_device_rad;
    double* g_device_z;
    double* g_device_tor;
    double* g_device_vrad;
    double* g_device_vz;
    double* g_device_vtor;
    int* g_device_detcellid;
    cudaMalloc((void **) &(g_device_rad), dimX);
    cudaMalloc((void **) &(g_device_z),   dimX);
    cudaMalloc((void **) &(g_device_tor), dimX);
    cudaMalloc((void **) &(g_device_vrad), dimX);
    cudaMalloc((void **) &(g_device_vz),   dimX);
    cudaMalloc((void **) &(g_device_vtor), dimX);
    cudaMalloc((void **) &(g_device_detcellid), dimXI);
    
    if (!FASTMODE){
        g_host->rad = (double*)malloc(dimX);
        g_host->z   = (double*)malloc(dimX);
        g_host->tor = (double*)malloc(dimX);
        g_host->vrad = (double*)malloc(dimX);
        g_host->vz   = (double*)malloc(dimX);
        g_host->vtor = (double*)malloc(dimX);
        g_host->detcellid = (int*)malloc(dimXI);
        load_beam(g_host, beam, shot, run);
        
        for (int i=0; i<run.block_size * run.block_number; ++i){
            g_host->detcellid[i] = -1;
        }    
        
        cudaMemcpy(&(g_device_rad),       &(g_host->vrad),      dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g_device_z),         &(g_host->z),         dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g_device_tor),       &(g_host->tor),       dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g_device_vrad),      &(g_host->vrad),      dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g_device_vz),        &(g_host->vz),        dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g_device_vtor),      &(g_host->vtor),      dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g_device_detcellid), &(g_host->detcellid), dimXI, cudaMemcpyHostToDevice);
    }else{
        BeamProfile dev_beam_prof;
        init_beam_profile(&dev_beam_prof, shot);
        printf("i84 \n");//# generate_coords <<< run.block_number, run.block_size >>> (*g, *s, beam, dev_beam_prof);
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

void coord_memcopy(TaigaGlobals *g_host, TaigaGlobals *g, TaigaCommons *s, BeamProp beam, ShotProp shot, RunProp run){
    size_t dimX = run.block_size * run.block_number * sizeof(double);
    size_t dimXI = run.block_size * run.block_number * sizeof(int);
    cudaMemcpy(&(g_host->rad),       &(g->rad),       dimX,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->z),         &(g->z),         dimX,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->tor),       &(g->tor),       dimX,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->vrad),      &(g->vrad),      dimX,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->vz),        &(g->vz),        dimX,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->vtor),      &(g->vtor),      dimX,  cudaMemcpyHostToDevice);
    cudaMemcpy(&(g_host->detcellid), &(g->detcellid), dimXI, cudaMemcpyHostToDevice);
}

void init_device_structs(TaigaGlobals *g_host, TaigaGlobals *g, TaigaCommons *s_host, TaigaCommons *s, BeamProp beam, ShotProp shot, RunProp run){
    g_host->particle_number = run.particle_number;
    s_host->max_step_number = run.step_device;        // N_step
    s_host->step_counter    = 0;
    s_host->eperm           = ELEMENTARY_CHARGE/ AMU/ beam.mass;
    s_host->timestep        = run.timestep;
    cudaMemcpy(&g, &g_host, sizeof(TaigaGlobals), cudaMemcpyHostToDevice);
    cudaMemcpy(&s, &s_host, sizeof(TaigaCommons), cudaMemcpyHostToDevice);
}

