#include <cuda.h>
#include "taiga_init.h"

void init_taiga(device_global *g, device_shared *s){
    size_t dimP = sizeof(double*);
    cudaMalloc((void **) &(s->spline_grid), 2*dimP);
    cudaMalloc((void **) &(s->bspline), 3*dimP);
    cudaMalloc((void **) &(s->espline), 3*dimP);
}

void init_device(device_global *g, device_shared *s){
    s->step_counter = 0;
    s->max_step_number = 0;
    s->eperm = 0;
    s-> timestep = 0;
    s->espline_on = false;
    g->particle_number = 0;
}

void init_grid(device_shared *s_host, device_shared *s, shot_prop shot){
    
    size_t dimGP = 2*sizeof(int*);
    cudaMalloc((void **) &s->grid_size, dimGP); //  cudaMalloc((void**) &, dim);  
    s_host->grid_size = (int*)malloc(dimGP);

    size_t dimG = 2*sizeof(double*);
    cudaMalloc((void **) &s->spline_grid, dimG);
    s_host->spline_grid = (double**)malloc(dimG);

    s_host->grid_size[0] = read_vector(&(s_host->spline_grid[0]), "input/fieldSpl", shot.name, "r.spline");
    size_t dimR = s_host->grid_size[0] * sizeof(double);
    cudaMalloc((void **) &(s->spline_grid[0]), dimR);
    
    s_host->grid_size[1] = read_vector(&(s_host->spline_grid[1]), "input/fieldSpl", shot.name, "z.spline");
    size_t dimZ = s_host->grid_size[1] * sizeof(double);
    cudaMalloc((void **) &(s->spline_grid[1]), dimZ);
    
//    size_t dimRZ = (dev_shared->grid_size[0]-1) * (dev_shared->grid_size[1]-1) * sizeof(double); // EZ grid

    cudaMemcpy(&(s->spline_grid[0]), &(s_host->spline_grid[0]), dimR, cudaMemcpyHostToDevice);
    cudaMemcpy(&(s->spline_grid[1]), &(s_host->spline_grid[1]), dimZ, cudaMemcpyHostToDevice);
    cudaMemcpy(&(s->spline_grid),    &(s_host->spline_grid),    dimG, cudaMemcpyHostToDevice);
}

void init_coords(device_global *g_host, device_global *g, beam_prop beam, shot_prop shot, run_prop run){

    size_t dimXP = 6*sizeof(double*);
    cudaMalloc((void**) &g_host->coords, dimXP);
    g_host->coords = (double**)malloc(dimXP);
    cudaMemcpy(&(g->coords), &(g_host->coords), dimXP, cudaMemcpyHostToDevice);

    size_t dimX = run.block_size * run.block_number * sizeof(double);

    for (int i=0; i<6; i++){
        cudaMalloc((void **) &g->coords[i], dimX);
    }

    if (!FASTMODE){
        //#load_beam(g_host.coords, beam, shot, run);

        for (int i=0; i<6; i++){
//            g_host->coords[i] = malloc(dimX);
            cudaMemcpy(&(g->coords[i]), &(g_host->coords[i]), dimX, cudaMemcpyHostToDevice);
        }
    }
}

void init_beam_profile(beam_profile *dev_prof, shot_prop shot){
    beam_profile host_prof;
    
//#    init_ion_profile(shot.name, host_prof);

    size_t dimBF = sizeof(beam_profile);
    size_t dimBR = sizeof(double)*host_prof.radial.N;
    size_t dimBX = sizeof(double)*host_prof.cross_section.N;

//#    load_ion_profile(shot.name, host_prof);

    cudaMalloc((void **) &(dev_prof->radial), dimBF);
    cudaMalloc((void **) &(dev_prof->radial.grid),    dimBR);
    cudaMalloc((void **) &(dev_prof->radial.profile), dimBR);
    cudaMemcpy(&(dev_prof->radial.grid),    &(host_prof.radial.grid),    dimBR, cudaMemcpyHostToDevice);
    cudaMemcpy(&(dev_prof->radial.profile), &(host_prof.radial.profile), dimBR, cudaMemcpyHostToDevice);

    cudaMalloc((void **) &(dev_prof->cross_section), dimBF);
    cudaMalloc((void **) &(dev_prof->cross_section.grid),    dimBX);
    cudaMalloc((void **) &(dev_prof->cross_section.profile), dimBX);
    cudaMemcpy(&(dev_prof->cross_section.grid),    &(host_prof.cross_section.grid),    dimBX, cudaMemcpyHostToDevice);
    cudaMemcpy(&(dev_prof->cross_section.profile), &(host_prof.cross_section.profile), dimBX, cudaMemcpyHostToDevice);    
}
