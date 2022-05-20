#include "cuda.h"
#include "utils/dataio/data_import.h"
#include "utils/taiga_constants.h"
#include "utils/physics.h"
#include "utils/basic_functions.c"
#include "utils/cuda.cu"

void init_host(TaigaGlobals *host_global, TaigaCommons *host_common){
    host_common->max_step_number = 0;
    host_common->eperm = 0;
    host_common->timestep = 0;
    host_common->ionisation_energy = INFINITY;
    host_common->is_electric_field_on = false;
    host_common->is_magnetic_field_perturbation = false;
    host_common->is_ionisation_on = false;
    host_global->particle_number = 0;
}

void init_grid(ShotProp shot, RunProp run, TaigaCommons *host_common, TaigaCommons *shared_common){
    int* shared_grid_size;
    double *shared_rgrid, *shared_zgrid;
    size_t size_commons = sizeof(TaigaCommons);
    size_t size_grid_dim = 2*sizeof(int);
    
    host_common->grid_size = (int*)malloc(size_grid_dim);
    if (run.field_interpolation_method == CUBIC_SPLINE) {
        host_common->grid_size[0] = read_vector(&host_common->spline_rgrid, "input/fieldSpl", shot.name, "r.spline");
        host_common->grid_size[1] = read_vector(&host_common->spline_zgrid, "input/fieldSpl", shot.name, "z.spline");
    }else if (run.field_interpolation_method == CUBIC_BSPLINE) {
        host_common->grid_size[0] = read_vector(&host_common->spline_rgrid, "input/fieldSpl", shot.name, "r.bspl");
        host_common->grid_size[1] = read_vector(&host_common->spline_zgrid, "input/fieldSpl", shot.name, "z.bspl");
    }else{
        printf("Error: invalid field interpolation method in init.\n");
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

    if (run.debug){
        printf(" GRID SIZE: %d %d \n", host_common->grid_size[0], host_common->grid_size[1]);
    }
}

void free_grid(TaigaCommons *shared_common){
    cudaFree(shared_common->grid_size);
    cudaFree(shared_common->spline_rgrid);
    cudaFree(shared_common->spline_zgrid);
}


void init_device_structs(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *shared_global, TaigaCommons *shared_common){
    shared_global->particle_number       = run.particle_number;
    shared_common->max_step_number       = run.step_device;
    shared_common->eperm                 = beam.charge * ELEMENTARY_CHARGE / AMU / get_mass(beam.species, beam.charge);
    shared_common->timestep              = run.timestep;
    shared_common->solver                = run.solver;
    shared_common->field_interpolation_method = run.field_interpolation_method;
    shared_common->is_ionisation_on      = run.is_ionisation_on;
    shared_common->ionisation_energy     = get_ionisation_energy(beam.species, beam.charge);
}

void set_particle_number(RunProp *run, TaigaGlobals *host_global, TaigaGlobals *shared_global){
    if (run->init_source == READ_COORDINATES){
        double *X_temp;
        host_global->particle_number = read_vector(&X_temp, "input", "manual_profile", "rad.dat");
        run->block_number = host_global->particle_number / run->block_size+1;
        free(X_temp);
    }else{
        host_global->particle_number = run->block_size * run->block_number;
    }
    shared_global->particle_number = host_global->particle_number;
}