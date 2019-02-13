inline void cErrorCheck(const char *file, int line) {
cudaThreadSynchronize();
cudaError_t err = cudaGetLastError();
if (err != cudaSuccess) {
    printf("Error: %s\n", cudaGetErrorString(err));
    printf(" @ %s: %d\n", file, line);
    exit(-1);
}
}

void set_cuda(int debug_flag){
    int num_devices, device_i, active_device=0;
    cudaGetDeviceCount(&num_devices);
    if (debug_flag) printf("Number of devices: %d\n", num_devices);
    
    if (num_devices > 1) {
        int max_multiprocessors = 0, active_device = 0;
        for (device_i = 0; device_i < num_devices; device_i++) {
            cudaDeviceProp properties;
            cudaGetDeviceProperties(&properties, device_i);
            if (max_multiprocessors < properties.multiProcessorCount) {
                max_multiprocessors = properties.multiProcessorCount;
                active_device = device_i;
            }
            if (debug_flag){
                printf("%d:%s\n",device_i,&properties.name);
                printf("\tL2Cache:\t%d", properties.l2CacheSize);
                printf("\tNumber of cores:\t%d", properties.warpSize);        
                printf("\tKernels:\t%d", properties.concurrentKernels);
                printf("\tThreads:\t%d", properties.maxThreadsPerMultiProcessor);
                printf("\tClock:\t%d", properties.clockRate/1024);
                printf("\n");
            }
        }
        cudaSetDevice(active_device);
        if (debug_flag){
            for (device_i = 0; device_i < num_devices; device_i++) {
                if(device_i==active_device) printf("-->");
                cudaDeviceProp properties;
                cudaGetDeviceProperties(&properties, device_i);
                printf("\t%d:\t%s\n",device_i,&properties.name);
            }
        }
    }

    cudaDeviceProp prop;
    cudaGetDevice(&active_device);
    cudaGetDeviceProperties(&prop, active_device);    
    printf("Active card:\t%s\n", &prop.name);
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
