#include <cuda.h>
#include "taiga_init.h"
#include "running/generate_coords.cuh"
#include "dataio/beam.h"

void init_taiga(taiga_globals *g, taiga_commons *s){
}

void init_device(taiga_globals *g, taiga_commons *s){
    s->step_counter = 0;
    s->max_step_number = 0;
    s->eperm = 0;
    s-> timestep = 0;
    s->espline_on = false;
    g->particle_number = 0;
}

void init_grid(shot_prop shot, run_prop run, taiga_commons *s_host, taiga_commons *s){
    
    double *host_rgrid, *dev_rgrid, *host_zgrid, *dev_zgrid;
    size_t dimGP = 2*sizeof(int);
    cudaMalloc((void **) &s->grid_size, dimGP); //  cudaMalloc((void**) &, dim);
    s_host->grid_size = (int*)malloc(dimGP);
    
    s_host->grid_size[0] = read_vector(&host_rgrid, "input/fieldSpl", shot.name, "r.spline");
    size_t dimR = s_host->grid_size[0] * sizeof(double);
    cudaMalloc((void **) &dev_rgrid, dimR);
    cudaMemcpy(dev_rgrid, host_rgrid, dimR, cudaMemcpyHostToDevice);
    s->spline_rgrid = dev_rgrid;

    s_host->grid_size[1] = read_vector(&host_zgrid, "input/fieldSpl", shot.name, "z.spline");
    size_t dimZ = s_host->grid_size[1] * sizeof(double);
    cudaMalloc((void **) &dev_zgrid, dimZ);
    cudaMemcpy(dev_zgrid, host_zgrid, dimZ, cudaMemcpyHostToDevice);
    s->spline_zgrid = dev_zgrid;
    
    cudaMemcpy(&(s->grid_size),    &(s_host->grid_size),    dimGP, cudaMemcpyHostToDevice);
    
    printf(" GRID SIZE: %d %d \n", s_host->grid_size[0], s_host->grid_size[1]);
}

void init_coords(taiga_globals *g_host, taiga_globals *g, taiga_commons *s, beam_prop beam, shot_prop shot, run_prop run){
    
    size_t dimX = run.block_size * run.block_number * sizeof(double);
    size_t dimXI = run.block_size * run.block_number * sizeof(int);
    
    cudaMalloc((void **) &(g->rad), dimX);
    cudaMalloc((void **) &(g->z),   dimX);
    cudaMalloc((void **) &(g->tor), dimX);
    cudaMalloc((void **) &(g->vrad), dimX);
    cudaMalloc((void **) &(g->vz),   dimX);
    cudaMalloc((void **) &(g->vtor), dimX);
    cudaMalloc((void **) &(g->detcellid), dimXI);
    
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
        
        cudaMemcpy(&(g->rad),       &(g_host->vrad),      dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g->z),         &(g_host->z),         dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g->tor),       &(g_host->tor),       dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g->vrad),      &(g_host->vrad),      dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g->vz),        &(g_host->vz),        dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g->vtor),      &(g_host->vtor),      dimX,  cudaMemcpyHostToDevice);
        cudaMemcpy(&(g->detcellid), &(g_host->detcellid), dimXI, cudaMemcpyHostToDevice);
    }else{
        beam_profile dev_beam_prof;
        init_beam_profile(&dev_beam_prof, shot);
        printf("i73 \n");//# generate_coords <<< run.block_number, run.block_size >>> (*g, *s, beam, dev_beam_prof);
    }
}

void init_beam_profile(beam_profile *dev_prof, shot_prop shot){
    /*#beam_profile host_prof;
    
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

void coord_memcopy(taiga_globals *g_host, taiga_globals *g, taiga_commons *s, beam_prop beam, shot_prop shot, run_prop run){
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

void init_device_structs(taiga_globals *g_host, taiga_globals *g, taiga_commons *s_host, taiga_commons *s, beam_prop beam, shot_prop shot, run_prop run){
    g_host->particle_number = run.particle_number;
    s_host->max_step_number = run.step_device;        // N_step
    s_host->step_counter    = 0;
    s_host->eperm           = ELEMENTARY_CHARGE/ AMU/ beam.mass;
    s_host->timestep        = run.timestep;
    cudaMemcpy(&g, &g_host, sizeof(taiga_globals), cudaMemcpyHostToDevice);
    cudaMemcpy(&s, &s_host, sizeof(taiga_commons), cudaMemcpyHostToDevice);
}

