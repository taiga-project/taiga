#include "interface/compiler_definition.h"
#include "utils/taiga_constants.h"

/*
DEFAULT_FLAGS = $(VERSION) -D'RENATE=0' -D'READINPUTPROF=1' -D'FASTMODE=0'
RENATE_FLAGS = $(VERSION) -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=0'
RENATE_FAST_FLAGS = $(VERSION) -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=1'
 */

void read_compiler_definition(RunProp *run) {
    if (RENATE == 1) {
        run->init_source = READ_RENATE_OD;
    } else {
        run->init_source = READ_COORDINATES;
    }

    switch(FASTMODE) {
        case 1:
            run->mode = JUST_WRITE; break;
        case 2:
            run->mode = MINIMAL_IO; break;
        default:
            run->mode = ALL_IO; break;
    }
}