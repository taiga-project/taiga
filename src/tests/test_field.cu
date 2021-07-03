#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#include "utils/taiga_constants.h"
#include "utils/prop.h"
#include "utils/basic_functions.h"
#include "dataio/data_import.c"
#include "dataio/field_import.cu"
#include "dataio/parameter_reader.c"
#include "init/sync.cu"
#include "init/init.cu"
#include "core/cyl2tor.cu"
#include "core/detection.cu"
#include "core/localise_field.cu"
#include "core/bspline.cu"


#define GRID_RES 33

__device__ double (*calculate_local_field)(TaigaCommons *c, const int *local_spline_indices,
                                           const double *local_spline, double dr, double dz);

__device__ double (*get_dr)(TaigaCommons *c, const int *local_spline_indices, double R);
__device__ double (*get_dz)(TaigaCommons *c, const int *local_spline_indices, double Z);

__global__ void calculate_field_grid(TaigaCommons *c, double *R, double *Z, double *field, double *psi_n){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    double r = R[idx];
    double z = Z[idx];
    
    int local_spline_indices[2];
    local_spline_indices[0] = SPLINE_INDEX_ERROR;
    local_spline_indices[1] = SPLINE_INDEX_ERROR;

    double local_spline_brad[16];
    double local_spline_bz[16];
    double local_spline_btor[16];

    double local_spline_erad[16];
    double local_spline_ez[16];
    double local_spline_etor[16];

    double local_spline_polflux[16];

    double local_brad=0, local_bz=0, local_btor=0;
    double local_erad=0, local_ez=0, local_etor=0;
    double local_psi_n = 0;
    double dr, dz;

    if (c->field_interpolation_method == CUBIC_SPLINE){
        get_coefficients = &get_coefficients_with_splines;
        calculate_local_field = &calculate_local_field_with_splines;
        get_dr = &get_dr_with_splines;
        get_dz = &get_dz_with_splines;
    }else if (c->field_interpolation_method == CUBIC_BSPLINE){
        get_coefficients = &get_coefficients_with_bsplines;
        calculate_local_field = &calculate_local_field_with_bsplines;
        get_dr = &get_dr_with_bsplines;
        get_dz = &get_dz_with_bsplines;
    }else{
            return;
    }

    copy_local_field(c, r, z, local_spline_indices,
                     local_spline_brad, local_spline_bz, local_spline_btor,
                     local_spline_erad, local_spline_ez, local_spline_etor,
                     local_spline_polflux);
    dr = (*get_dr)(c, local_spline_indices, r);
    dz = (*get_dz)(c, local_spline_indices, z);
    local_brad = (*calculate_local_field)(c, local_spline_indices, local_spline_brad, dr, dz);
    local_bz   = (*calculate_local_field)(c, local_spline_indices, local_spline_bz,   dr, dz);
    local_btor = (*calculate_local_field)(c, local_spline_indices, local_spline_btor, dr, dz);
    local_psi_n = (*calculate_local_field)(c, local_spline_indices, local_spline_polflux, dr, dz);
    field[idx] = local_brad;
    field[idx+GRID_RES*GRID_RES] = local_bz;
    field[idx+2*GRID_RES*GRID_RES] = local_btor;
    psi_n[idx] = local_psi_n;
}

void test_field(int field_interpolation_method){
    char* field_interpolation_name;
    switch (field_interpolation_method) {
        case CUBIC_SPLINE:
            field_interpolation_name = "spline";
            break;
        case CUBIC_BSPLINE:
            field_interpolation_name = "bspline";
            break;
        default:
            printf("Invalid interpolation value\n");
            exit(1);
    }

    ShotProp shot; init_shot_prop(&shot);
    BeamProp beam; init_beam_prop(&beam);
    RunProp run;   init_run_prop(&run);

    set_cuda(run.debug);
    
    parameter_reader(&beam, &shot, &run);
    runnumber_reader(&shot, &run);
    
    TaigaGlobals *device_global, *host_global, *shared_global;
    TaigaCommons *device_common, *host_common, *shared_common;
    
    size_t size_global = sizeof(TaigaGlobals);
    size_t size_commons = sizeof(TaigaCommons);
    
    host_global = (TaigaGlobals*)malloc(size_global);
    shared_global = (TaigaGlobals*)malloc(size_global);
    host_common = (TaigaCommons*)malloc(size_commons);
    shared_common = (TaigaCommons*)malloc(size_commons);
    
    cudaMalloc((void **) &device_global, size_global);
    cudaMalloc((void **) &device_common, size_commons);

    init_host(host_global, host_common);
    run.field_interpolation_method = field_interpolation_method;
    run.is_magnetic_field_perturbation = true;

    init_grid(shot, run, host_common, shared_common);

    switch (field_interpolation_method) {
        case CUBIC_SPLINE:
            magnetic_field_read_and_init(shot, run, host_common, shared_common);
            poloidal_flux_read_and_init(shot, run, host_common, shared_common);
            break;
        case CUBIC_BSPLINE:
            magnetic_field_read_and_init_with_bsplines(shot, run, host_common, shared_common);
            poloidal_flux_read_and_init_with_bsplines(shot, run, host_common, shared_common);
            break;
        default:
            printf("Invalid interpolation value\n");
            exit(1);
    }


    init_device_structs(beam, shot, run, shared_global, shared_common);
    sync_device_structs(device_global, shared_global, device_common, shared_common);

    double *host_field, *device_field;
    double *host_R, *device_R;
    double *host_Z, *device_Z;
    double *host_psi_n, *device_psi_n;
    long grid_size = GRID_RES*GRID_RES;
    size_t dim_tmp = sizeof(double)*grid_size;
    host_field = (double *) malloc(3*dim_tmp);
    host_R = (double *) malloc(dim_tmp);
    host_Z = (double *) malloc(dim_tmp);
    host_psi_n = (double *) malloc(dim_tmp);
    cudaMalloc((void **) &(device_field), 3*dim_tmp);
    cudaMalloc((void **) &(device_R), dim_tmp);
    cudaMalloc((void **) &(device_Z), dim_tmp);
    cudaMalloc((void **) &(device_psi_n), dim_tmp);

    double R_max=0.8;
    double R_min=0.3;
    double Z_max=0.4;
    double Z_min=-0.4;

    for(int i=0; i<GRID_RES; ++i){
        for(int j=0; j<GRID_RES; ++j){
            int index=i*GRID_RES+j;
            host_R[index] = i*(R_max-R_min)/(GRID_RES-1)+R_min;
            host_Z[index] = j*(Z_max-Z_min)/(GRID_RES-1)+Z_min;
            host_field[index]=UNDEFINED_FLOAT;
            host_field[grid_size+index]=UNDEFINED_FLOAT;
            host_field[2*grid_size+index]=UNDEFINED_FLOAT;
            host_psi_n[index]=UNDEFINED_FLOAT;
        }
    }

    cudaMemcpy(device_field, host_field, 3*dim_tmp, cudaMemcpyHostToDevice);
    cudaMemcpy(device_R, host_R, dim_tmp, cudaMemcpyHostToDevice);
    cudaMemcpy(device_Z, host_Z, dim_tmp, cudaMemcpyHostToDevice);
    cudaMemcpy(device_psi_n, host_psi_n, dim_tmp, cudaMemcpyHostToDevice);

    calculate_field_grid <<< GRID_RES, GRID_RES >>> (device_common, device_R, device_Z, device_field, device_psi_n);

    cudaMemcpy(host_field, device_field, 3*dim_tmp, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_psi_n, device_psi_n, dim_tmp, cudaMemcpyDeviceToHost);

    FILE *fp;
    fp = fopen (concat(run.folder_out, "/test_", field_interpolation_name, "_brad.dat", NULL), "w");
    for (int i=0; i<grid_size; ++i){
        fprintf(fp, "%lf %lf %.18lg\n", host_R[i], host_Z[i], host_field[i]);
    }
    fclose(fp);
    
    fp = fopen (concat(run.folder_out, "/test_", field_interpolation_name, "_bz.dat", NULL), "w");
    for (int i=0; i<grid_size; ++i){
        fprintf(fp, "%lf %lf %.18lg\n", host_R[i], host_Z[i], host_field[grid_size+i]);
    }
    fclose(fp);
    
    fp = fopen (concat(run.folder_out, "/test_", field_interpolation_name, "_btor.dat", NULL), "w");
    for (int i=0; i<grid_size; ++i){
        fprintf(fp, "%lf %lf %.18lg\n", host_R[i], host_Z[i], host_field[2*grid_size+i]);
    }
    fclose(fp);

    fp = fopen (concat(run.folder_out, "/test_", field_interpolation_name, "_psi.dat", NULL), "w");
    for (int i=0; i<grid_size; ++i){
        fprintf(fp, "%lf %lf %.18lg\n", host_R[i], host_Z[i], host_psi_n[i]);
    }
    fclose(fp);
    printf("Data exported to: %s/test_%s_*.dat\n", run.folder_out, field_interpolation_name);
    cudaDeviceReset() ;
}

int main(){
    test_field(CUBIC_SPLINE);
    test_field(CUBIC_BSPLINE);
}