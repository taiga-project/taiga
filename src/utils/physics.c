#include "taiga_constants.h"
#include "physics.h"

//https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
double get_mass(char *name_of_ion, int charge){
    if (!strcmp(name_of_ion,"D") || !strcmp(name_of_ion,"H2")){
        return 2.0141017781212 - charge * MASS_OF_ELECTRON;
    }else if (!strcmp(name_of_ion,"Li") || !strcmp(name_of_ion,"Li7")){
        return 7.016003436645 - charge * MASS_OF_ELECTRON;
    }else if (!strcmp(name_of_ion,"Na") || !strcmp(name_of_ion,"Na23")){
        return 22.989769282019 - charge * MASS_OF_ELECTRON;
    }else if (!strcmp(name_of_ion,"K") || !strcmp(name_of_ion,"K39")){
        return 38.963706486449 - charge * MASS_OF_ELECTRON;
    }else if (!strcmp(name_of_ion,"Rb") || !strcmp(name_of_ion,"Rb85")){
        return 84.911789737954 - charge * MASS_OF_ELECTRON;
    }else if (!strcmp(name_of_ion,"Cs") || !strcmp(name_of_ion,"Cs133")){
        return 132.905451961080 - charge * MASS_OF_ELECTRON;
    }else{
        try{
            return atof(name_of_ion);
        }catch (...){
            printf("Error: Invalid ion species: %s\n", name_of_ion);
            exit(1);
        }
    }
}

//http://www-personal.umich.edu/~cowley/ionen.htm
double get_ionisation_energy(char *name_of_ion, int charge){
    switch (charge) {
        case 0:
            if (!strcmp(name_of_ion, "Li") || !strcmp(name_of_ion, "Li7")) {
                return 5.39171;
            } else if (!strcmp(name_of_ion, "Na") || !strcmp(name_of_ion, "Na23")) {
                return 5.13908;
            } else if (!strcmp(name_of_ion, "K") || !strcmp(name_of_ion, "K39")) {
                return 4.34066;
            } else if (!strcmp(name_of_ion, "Rb") || !strcmp(name_of_ion, "Rb85")) {
                return 4.17713;
            } else if (!strcmp(name_of_ion, "Cs") || !strcmp(name_of_ion, "Cs133")) {
                return 3.89390;
            } else {
                printf("Warning: Invalid ion species: %s\n         Secondary ionisation module turned off", name_of_ion);
                return INFINITY;
            }
        case 1:
            if (!strcmp(name_of_ion, "Li") || !strcmp(name_of_ion, "Li7")) {
                return 75.64009;
            } else if (!strcmp(name_of_ion, "Na") || !strcmp(name_of_ion, "Na23")) {
                return 47.2864;
            } else if (!strcmp(name_of_ion, "K") || !strcmp(name_of_ion, "K39")) {
                return 31.63;
            } else if (!strcmp(name_of_ion, "Rb") || !strcmp(name_of_ion, "Rb85")) {
                return 27.2895;
            } else if (!strcmp(name_of_ion, "Cs") || !strcmp(name_of_ion, "Cs133")) {
                return 23.15744;
            } else {
                printf("Warning: Invalid ion species: %s\n         Secondary ionisation module turned off", name_of_ion);
                return INFINITY;
            }
        default:
            printf("Warning: Invalid charge number: %d\n         Secondary ionisation module turned off", charge);
            return INFINITY;
    }
}
