#include <cuda.h>
#include "utils/taiga_constants.h"
#include "utils/basic_functions.c"

void init_host(TaigaGlobals *host_global, TaigaCommons *host_common){
    host_common->step_counter = 0;
    host_common->max_step_number = 0;
    host_common->eperm = 0;
    host_common->timestep = 0;
    host_common->is_electric_field_on = false;
    host_common->is_magnetic_field_perturbation = false;
    host_global->particle_number = 0;
}

void init_grid(ShotProp shot, RunProp run, TaigaCommons *host_common, TaigaCommons *shared_common){
    int* shared_grid_size;
    double *shared_rgrid, *shared_zgrid;
    size_t size_commons = sizeof(TaigaCommons);
    size_t size_grid_dim = 2*sizeof(int);
    
    host_common->grid_size = (int*)malloc(size_grid_dim);
    if (host_common->field_interpolation_method == CUBIC_SPLINE) {
        host_common->grid_size[0] = read_vector(&host_common->spline_rgrid, "input/fieldSpl", shot.name, "r.spline");
        host_common->grid_size[1] = read_vector(&host_common->spline_zgrid, "input/fieldSpl", shot.name, "z.spline");
    }else if (host_common->field_interpolation_method == CUBIC_BSPLINE) {
        host_common->grid_size[0] = read_vector(&host_common->spline_rgrid, "input/fieldSpl", shot.name, "r.bspline");
        host_common->grid_size[1] = read_vector(&host_common->spline_zgrid, "input/fieldSpl", shot.name, "z.bspline");
    }else{
        printf("Error: invalid field interpolation method in init.\n");
        exit(1);
    }
    size_t size_R = host_common->grid_size[0] * sizeof(double);
    size_t size_Z = host_common->grid_size[1] * sizeof(double);
    
    memcpy(shared_common, host_common, size_commons);
    
    cudaMalloc((void **) &shared_grid_size, size_grid_dim);
    cudaMemcpy(shared_grid_size, host_common->grid_size, size_grid_dim, cudaMemcpyHostToDevice);
    shared_common->grid_size = shared_grid_size;
    
    cudaMalloc((void **) &shared_rgrid, size_R);
    cudaMemcpy(shared_rgrid, host_common->spline_rgrid, size_R, cudaMemcpyHostToDevice);
    shared_common->spline_rgrid = shared_rgrid;
    
    cudaMalloc((void **) &shared_zgrid, size_Z);
    cudaMemcpy(shared_zgrid, host_common->spline_zgrid, size_Z, cudaMemcpyHostToDevice);
    shared_common->spline_zgrid = shared_zgrid;
    
    printf(" GRID SIZE: %d %d \n", host_common->grid_size[0], host_common->grid_size[1]);
}

void init_device_structs(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *shared_global, TaigaCommons *shared_common){
    shared_global->particle_number       = run.particle_number;
    shared_common->max_step_number       = run.step_device;
    shared_common->step_counter          = 0;
    shared_common->eperm                 = ELEMENTARY_CHARGE / AMU / beam.mass;
    shared_common->timestep              = run.timestep;
    shared_common->solver                = run.solver;
    shared_common->field_interpolation_method = run.field_interpolation_method;
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