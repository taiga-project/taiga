#include "test_basic_functions.c"
#include "test_taiga_init.cu"

#define LENGTH_TMP 10

int main(){
    test_concat();
    
    test_init_grid();
    test_init_coords();
    test_init_detector();
    test_init_detector_full();
    test_detector_conversion();
}

