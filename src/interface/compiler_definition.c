#include "interface/compiler_definition.h"
#include "utils/taiga_constants.h"

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