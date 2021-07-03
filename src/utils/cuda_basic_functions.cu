#include "cuda_basic_functions.cuh"

inline void cErrorCheck(const char *file, int line){
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess){
        printf("Error: %s\n", cudaGetErrorString(err));
        printf(" @ %s: %d\n", file, line);
        exit(1);
    }
}

void set_cuda(int debug_flag){
    int num_devices, device_i, active_device=0;
    cudaDeviceProp properties;
    cudaGetDeviceCount(&num_devices);
    if (debug_flag) printf("Number of devices: %d\n", num_devices);

    if (num_devices > 1 || debug_flag) {
        int max_multiprocessors = 0;
        for (device_i = 0; device_i < num_devices; ++device_i){
            cudaGetDeviceProperties(&properties, device_i);
            if (max_multiprocessors < properties.multiProcessorCount){
                max_multiprocessors = properties.multiProcessorCount;
                active_device = device_i;
            }
        }
        cudaSetDevice(active_device);
        if (debug_flag){
            for (device_i = 0; device_i < num_devices; ++device_i){
                cudaGetDeviceProperties(&properties, device_i);
                if(device_i==active_device) printf("[*] "); else    printf("[ ] ");
                printf("Card %d:    %s\n", device_i, &properties.name);
                printf("      L2Cache:         %d\n", properties.l2CacheSize);
                printf("      Number of cores: %d\n", properties.warpSize);
                printf("      Kernels:         %d\n", properties.concurrentKernels);
                printf("      Threads:         %d\n", properties.maxThreadsPerMultiProcessor);
                printf("      Clock:           %d\n", properties.clockRate/1024);
                printf("\n");
            }
        }
    }

    cudaGetDevice(&active_device);
    cudaGetDeviceProperties(&properties, active_device);
    printf("Active card:\t%s\n", &properties.name);
}