#include <cuda.h>

void sync_device_structs(TaigaGlobals *g_device, TaigaGlobals *g_shared, TaigaCommons *c_device, TaigaCommons *c_shared){
    cudaMemcpy(c_device, c_shared, sizeof(TaigaCommons), cudaMemcpyHostToDevice);
    if (!FASTMODE){
        cudaMemcpy(g_device, g_shared, sizeof(TaigaGlobals), cudaMemcpyHostToDevice);
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
    int* host_detcellid =(int*)malloc(size_detcellid);

    cudaMemcpy(host_rad,  g_shared->rad,  size_coord, cudaMemcpyDeviceToHost); g_host->rad = host_rad;
    cudaMemcpy(host_z,    g_shared->z,    size_coord, cudaMemcpyDeviceToHost); g_host->z = host_z;
    cudaMemcpy(host_tor,  g_shared->tor,  size_coord, cudaMemcpyDeviceToHost); g_host->tor = host_tor;
    cudaMemcpy(host_vrad, g_shared->vrad, size_coord, cudaMemcpyDeviceToHost); g_host->vrad = host_vrad;
    cudaMemcpy(host_vz,   g_shared->vz,   size_coord, cudaMemcpyDeviceToHost); g_host->vz = host_vz;
    cudaMemcpy(host_vtor, g_shared->vtor, size_coord, cudaMemcpyDeviceToHost); g_host->vtor = host_vtor;
    cudaMemcpy(host_detcellid, g_shared->detcellid, size_detcellid, cudaMemcpyDeviceToHost); g_host->detcellid = host_detcellid;
}