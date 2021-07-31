#ifndef TEST_TAIGA_INIT_CUH
#define TEST_TAIGA_INIT_CUH

__global__ void test_init_grid_cuda(TaigaCommons* c, double *tmp);
__global__ void test_init_coords_cuda(TaigaGlobals* g, double *tmp);
void test_init_grid();
void test_init_coords();
__global__ void test_detector_cuda(DetectorProp *d,  double *tmp);
void test_init_detector();
__global__ void test_detector_struct(TaigaCommons *c, DetectorProp *d, double *tmp);
__global__ void test_detector_detcellid(TaigaGlobals *g, double *tmp);
__global__ void test_detector_full_cuda(DetectorProp *d, double *tmp);
void test_init_detector_full();
__global__ void test_detector_index_cuda(double *tmp);
void test_detector_conversion();
__global__ void test_renate_fast_cuda(TaigaGlobals* g, double *tmp);
void test_renate_fast(TaigaGlobals *g);

#endif //TEST_TAIGA_INIT_CUH
