#include "cuda.h"
#include "utils/cuda.cuh"

void sync_device_structs(TaigaGlobals *g_device, TaigaGlobals *g_shared, TaigaCommons *c_device, TaigaCommons *c_shared,
                         bool is_all_io){
    CHECK_ERROR(cudaMemcpy(c_device, c_shared, sizeof(TaigaCommons), cudaMemcpyHostToDevice));
    if (is_all_io){
        CHECK_ERROR(cudaMemcpy(g_device, g_shared, sizeof(TaigaGlobals), cudaMemcpyHostToDevice));
    }
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
    double* host_intensity =(double*)malloc(size_coord);
    double* host_time_of_flight =(double*)malloc(size_coord);
    int* host_detcellid =(int*)malloc(size_detcellid);

    CHECK_ERROR(cudaMemcpy(host_rad,  g_shared->rad,  size_coord, cudaMemcpyDeviceToHost));
    g_host->rad = host_rad;
    CHECK_ERROR(cudaMemcpy(host_z,    g_shared->z,    size_coord, cudaMemcpyDeviceToHost));
    g_host->z = host_z;
    CHECK_ERROR(cudaMemcpy(host_tor,  g_shared->tor,  size_coord, cudaMemcpyDeviceToHost));
    g_host->tor = host_tor;
    CHECK_ERROR(cudaMemcpy(host_vrad, g_shared->vrad, size_coord, cudaMemcpyDeviceToHost));
    g_host->vrad = host_vrad;
    CHECK_ERROR(cudaMemcpy(host_vz,   g_shared->vz,   size_coord, cudaMemcpyDeviceToHost));
    g_host->vz = host_vz;
    CHECK_ERROR(cudaMemcpy(host_vtor, g_shared->vtor, size_coord, cudaMemcpyDeviceToHost));
    g_host->vtor = host_vtor;
    CHECK_ERROR(cudaMemcpy(host_intensity, g_shared->intensity, size_coord, cudaMemcpyDeviceToHost));
    g_host->intensity = host_intensity;
    CHECK_ERROR(cudaMemcpy(host_time_of_flight, g_shared->time_of_flight, size_coord, cudaMemcpyDeviceToHost));
    g_host->time_of_flight = host_time_of_flight;
    CHECK_ERROR(cudaMemcpy(host_detcellid, g_shared->detcellid, size_detcellid, cudaMemcpyDeviceToHost));
    g_host->detcellid = host_detcellid;
}