#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parameter_reader.h"

#define MAXCHAR 1000

char* clean_string (char* str_in){
    char* str_out = str_in;
    ++str_out;    
    str_out[strlen(str_out)-2] = 0;
    return str_out;
}

void init_taiga_props(char* par_name, char* par_value, BeamProp *beam, ShotProp *shot, RunProp *run){
    
    double par_value_lf;    sscanf(par_value, "%lf", &par_value_lf);
    int par_value_d;        sscanf(par_value, "%d", &par_value_d);
    
    if (!strcmp(par_name, "shotnumber"))                    strcpy(shot->shotnumber, clean_string(par_value));
    else if (!strcmp(par_name, "time"))                     strcpy(shot->time, clean_string(par_value));
    else if (!strcmp(par_name, "runnumber"))                strcpy(run->runnumber, clean_string(par_value));
    else if (!strcmp(par_name, "matter") || !strcmp(par_name, "beammatter") || !strcmp(par_name, "species")){
        strcpy(beam->matter, clean_string(par_value));
        beam->mass = get_mass(beam->matter);
    }
    else if (!strcmp(par_name, "energy"))                   beam->energy = par_value_lf;
    else if (!strcmp(par_name, "vertical_deflection"))      beam->vertical_deflection = par_value_lf/180.0*PI;
    else if (!strcmp(par_name, "toroidal_deflection"))      beam->toroidal_deflection = par_value_lf/180.0*PI;
    else if (!strcmp(par_name, "angle")){                   // duplication (former name until 1.3.0)
        beam->toroidal_deflection = par_value_lf/180.0*PI;
        printf("Warning: Please rename 'angle' variable name to 'toroidal_deflection' in 'parameters.sh'!\n");
    }
    else if (!strcmp(par_name, "diameter"))                 beam->diameter = par_value_lf/1000.0;
    else if (!strcmp(par_name, "detector"))                 strcpy(shot->detector_geometry, clean_string(par_value));
    else if (!strcmp(par_name, "detector_mask"))            strcpy(shot->detector_mask, clean_string(par_value));
    else if (!strcmp(par_name, "electric_field_module"))    shot->is_electric_field_on = (bool)par_value_d;
    else if (!strcmp(par_name, "timestep"))                 run->timestep = par_value_lf;
    else if (!strcmp(par_name, "step_device"))              run->step_device = par_value_d;
    else if (!strcmp(par_name, "step_host"))                run->step_host = par_value_d;
    else if (!strcmp(par_name, "particles")){
        run->particle_number = (long)par_value_d;
        run->block_number    = par_value_d/run->block_size+1;
    }
}

int parameter_reader(BeamProp *beam, ShotProp *shot, RunProp *run){
    
    FILE *fp;
    char str[MAXCHAR];
    char* par_name;
    char* par_value;
    
    fp = fopen(run->parameter_file, "r");
    if (fp == NULL){
        printf("Could not open file %s",run->parameter_file);
        return 1;
    }
    while (fgets(str, MAXCHAR, fp) != NULL){
        if (str[0] != '#'){        
            par_name = strtok(str, "=");
            par_value = strtok(NULL, "#");
            init_taiga_props(par_name, par_value, beam, shot, run);
        }
    }
    fclose(fp);
    strcpy(shot->name, concat(shot->shotnumber, "_", shot->time, NULL));
    strcpy(run->folder_out, concat("results/", shot->shotnumber, "_", shot->time, NULL));
    return 0;
}

int runnumber_reader(ShotProp *shot, RunProp *run){
    if (strcmp(run->runnumber,"0") == 0){
        FILE *fp;
        char txt[10], *runnumber;
        
        fp = fopen(run->runnumber_file, "r");
        if (fp == NULL){
            printf("Could not open file %s",run->runnumber_file);
            return 1;
        }
        fgets(txt, 10, fp);
        runnumber = strtok(txt, "\n");
        strcpy(run->runnumber, runnumber);
    }
}

double get_mass(char *s){
    double mass;
    
    if (strcmp(s,"D")==0){
        mass = 2.013553212724;
    }else if (strcmp(s,"Li")==0){
        mass = 7.016004558;
    }else if (strcmp(s,"Na")==0){
        mass = 22.98976928;
    }else if (strcmp(s,"K")==0){
        mass = 39.9639984821;
    }else if (strcmp(s,"H2")==0){
        mass = 2.013553212724;
    }else if (strcmp(s,"Li7")==0){
        mass = 7.016004558;
    }else if (strcmp(s,"Na23")==0){
        mass = 22.98976928;
    }else if (strcmp(s,"K40")==0){
        mass = 39.9639984821;
    }else{
        try{
            mass = atof(s);
        }catch (...){
            mass = 7.016004558;
        }
    }
    return mass;
}
