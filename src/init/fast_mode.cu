#include <cuda.h>
#include "core/generate_coords.cuh"

void init_fastmode(BeamProp beam, ShotProp shot, RunProp run, TaigaGlobals *device_global){
    BeamProfile *device_prof;
    size_t size_prof = sizeof(BeamProfile);
    cudaMalloc((void **) &device_prof, size_prof);
    init_beam_profile(device_prof, shot);
    generate_coords <<< run.block_number, run.block_size >>> (device_global, beam, device_prof);
}