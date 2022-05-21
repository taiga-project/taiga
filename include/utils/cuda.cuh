#ifndef TAIGA_CUDA_CUH
#define TAIGA_CUDA_CUH

#define CHECK_ERROR(cuda_code) { manage_cuda_error((cuda_code), __FILE__, __LINE__); }

void set_cuda(int debug_flag);

inline void manage_cuda_error(cudaError_t cuda_code, const char *filename, int line) {
    if (cuda_code != cudaSuccess) {
        fprintf(stderr,"CUDA ERROR in %s (line %d):\n %s\n", filename, line, cudaGetErrorString(cuda_code));
    }
}

#endif //TAIGA_CUDA_CUH
