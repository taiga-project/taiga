#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include <time.h>
#include <math.h>
#include <string.h>

int ConvertSMVer2Cores(int major, int minor)
{
    // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
    typedef struct {
            int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
            int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] = {
        { 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
        { 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
        { 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
        { 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
        { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
        { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
        { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
        { 0x32, 192}, // Kepler Generation (SM 3.2) GK10x class
        { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
        { 0x37, 192}, // Kepler Generation (SM 3.7) GK21x class
        { 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
        { 0x52, 128}, // Maxwell Generation (SM 5.2) GM20x class
        {   -1, -1 }
    };

    int index = 0;
    while (nGpuArchCoresPerSM[index].SM != -1) {
            if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor) ) {
                    return nGpuArchCoresPerSM[index].Cores;
            }
            ++index;
    }
    printf("MapSMtoCores SM %d.%d is undefined (please update to the latest SDK)!\n", major, minor);
    return -1;
}

void DisplayHeader()
{
    const int kb = 1024;
    const int mb = kb * kb;

    int devCount;
    cudaGetDeviceCount(&devCount);
    //wcout << "CUDA Devices: " << endl << endl;
    printf("CUDA Devices: totally %d devices\n\n",devCount);
    for(int i = 0; i < devCount && i<20; ++i)
    {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        //printf("%d: %s\n",i,&props.name);
     //   wcout << i << ": " << props.name << ": " << props.major << "." << props.minor << endl;
         printf("%d: %s: %d.%d\n",i,props.name,props.major,props.minor); 
        printf("  CUDA Cores:\t\t%d (%dx%d)\n",
               ConvertSMVer2Cores(props.major, props.minor) * props.multiProcessorCount,
               props.multiProcessorCount,
               ConvertSMVer2Cores(props.major, props.minor));

        printf("  Global memory:\t%lf MiB\n", (double)props.totalGlobalMem / mb );
        printf("  Shared memory:\t%lf kiB\n" , (double)props.sharedMemPerBlock / kb );
        printf("  Constant memory:\t%lf kiB\n" , (double)props.totalConstMem / kb );
        printf("  Block registers:\t%d\n" , props.regsPerBlock );
        printf("  Multiprocessors:\t%d\n", props.multiProcessorCount);
        printf("  Warp size:\t\t%d\n", props.warpSize);
        printf("  Threads per block:\t%d\n", props.maxThreadsPerBlock);
        printf("  Max block dimensions:\t[%d, %d, %d]\n", props.maxThreadsDim[0] , props.maxThreadsDim[1] , props.maxThreadsDim[2]);
        printf("  Max grid dimensions:\t[%d, %d, %d]\n", props.maxGridSize[0] ,  props.maxGridSize[1] , props.maxGridSize[2]);
        printf("  Timeout (1:enabled):\t%d", props.kernelExecTimeoutEnabled);
        printf("\n");
    }
}

int main(){
    DisplayHeader();
    return 0;
}


