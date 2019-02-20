#if defined(_WIN32)
    #include <windows.h>
#else
    #include <unistd.h>
    void CopyFile(char* source, char* target, int sw){
        execl("/bin/cp", "-p", source, target);
        printf("Copy from %s to %s", source, target);
    }
#endif

inline void cErrorCheck(const char *file, int line){
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess){
        printf("Error: %s\n", cudaGetErrorString(err));
        printf(" @ %s: %d\n", file, line);
        exit(-1);
    }
}

void set_cuda(int debug_flag){
    int num_devices, device_i, active_device=0;
    cudaDeviceProp properties;
    cudaGetDeviceCount(&num_devices);
    if (debug_flag) printf("Number of devices: %d\n", num_devices);
    
    if (num_devices > 1 || debug_flag) {
        int max_multiprocessors = 0;
        for (device_i = 0; device_i < num_devices; device_i++){
            cudaGetDeviceProperties(&properties, device_i);
            if (max_multiprocessors < properties.multiProcessorCount){
                max_multiprocessors = properties.multiProcessorCount;
                active_device = device_i;
            }
        }
        cudaSetDevice(active_device);
        if (debug_flag){
            for (device_i = 0; device_i < num_devices; device_i++){
                cudaGetDeviceProperties(&properties, device_i);                
                if(device_i==active_device) printf("[*] "); else    printf("[ ] ");
                printf("Card %d:    %s\n",device_i,&properties.name);
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

double linear_interpolate(double *x_vector, int x_length, double *y_vector, int y_length, double x_value){
    int i;
    if (x_length != y_length)   printf("ERROR: in interpolation. Two input vectors have different length.");    
    for (i=1; (i<x_length) && (x_vector[i-1]>x_value); i++);    
    if(i>1){--i;}else{i=1;}    
    return y_vector[i] - (y_vector[i]-y_vector[i-1])*(x_value-x_vector[i-1])/(x_vector[i]-x_vector[i-1]);
}

int get_array_size(double *array){
    return (int)(sizeof(array)/sizeof(array[0]));
}

char* concat(const char *s1, const char *s2){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

char* concat(const char *s1, const char *s2, const char *s3){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+1);
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    return result;
}

char* concat(const char *s1, const char *s2, const char *s3, const char *s4){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+strlen(s4)+1);
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    strcat(result, s4);
    return result;
}

char* concat(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+strlen(s4)+strlen(s5)+1);
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    strcat(result, s4);
    strcat(result, s5);
    return result;
}

char* concat(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5, const char *s6){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+strlen(s4)+strlen(s5)+strlen(s6)+1);
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    strcat(result, s4);
    strcat(result, s5);
    strcat(result, s6);
    return result;
}

char* concat(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5, const char *s6, const char *s7){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+strlen(s4)+strlen(s5)+strlen(s6)+strlen(s7)+1);
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    strcat(result, s4);
    strcat(result, s5);
    strcat(result, s6);
    strcat(result, s7);
    return result;
}
