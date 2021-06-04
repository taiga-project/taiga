#include <cuda.h>
#include "utils/taiga_constants.h"
#include "utils/basic_functions.c"

void init_host(TaigaGlobals *g, TaigaCommons *s){
    s->step_counter = 0;
    s->max_step_number = 0;
    s->eperm = 0;
    s->timestep = 0;
    s->is_electric_field_on = false;
    s->is_magnetic_field_perturbation = false;
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