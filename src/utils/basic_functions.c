#include <stdarg.h>
#include <ctype.h>
#include "basic_functions.h"

#if defined(_WIN32)
    #include <windows.h>
    void mkdir(char *path, mode_t mode){
        _mkdir(path);
    }
#else
    #include <unistd.h>
    #include <sys/stat.h>
    #include <sys/types.h>
    void CopyFile(char* source, char* target, int sw){
        if(sw)  system(concat("cp -n", source, " ", target, NULL));
        else    system(concat("cp ", source, " ", target, NULL));
    }
#endif

void CopyFile(char* source, char* target){
    CopyFile(source, target, 0);
}

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

double linear_interpolate(double *x_vector, long x_length, double *y_vector, long y_length, double x_value){
    int i;
    if (x_length != y_length)   printf("ERROR: in interpolation. Two input vectors have different length.");
    for (i=1; (i<x_length) && (x_vector[i-1]>x_value); ++i);
    if(i>1){--i;}else{i=1;}
    return y_vector[i] - (y_vector[i]-y_vector[i-1])*(x_value-x_vector[i-1])/(x_vector[i]-x_vector[i-1]);
}

int get_array_size(double *array){
    return (int)(sizeof(array)/sizeof(array[0]));
}

void init_dir(char *folder, char *runnumber){ 
    init_dir(folder, runnumber, "");
}
void init_dir(char *folder, char *runnumber, char *subdir){
    mkdir(folder, 0777);
    mkdir(concat(folder, "/", runnumber, NULL), 0777);
    mkdir(concat(folder, "/", runnumber, "/", subdir, NULL), 0777);
}

char* concat(const char *s1, ...){
    const char *s;
    va_list args;
    char *r;//r[STRING_LENGTH];
    size_t arg_length = 1;
    
    va_start(args, s1);
    for(s=s1; s!=NULL; s=va_arg(args, const char*)){
        arg_length += strlen(s);
    }
    va_end(args);
    r = (char*)malloc(arg_length);
    r[0] = NULL;
    va_start(args, s1);
    for(s=s1; s!=NULL; s=va_arg(args, const char*)){
        strcat(r, s);
    }
    va_end(args);
    return r;
}

void string_to_lowercase(char* str){
    for(char *p=str; *p; ++p) *p=tolower(*p);
}