#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parameter_editor.h"

#define MAXCHAR 1000

char* clean_string (char* str_in){
    char* str_out = str_in;
    str_out++;    
    str_out[strlen(str_out)-2] = 0;
    return str_out;
}

void init_taiga_props(char* par_name, char* par_value, shot_prop *shot, beam_prop *beam, run_prop *run){
    
    double par_value_lf;    sscanf(par_value, "%lf", &par_value_lf);
    int par_value_d;        sscanf(par_value, "%d", &par_value_d);
        
    if (!strcmp(par_name, "shotnumber"))                    strcpy(shot->shotnumber, clean_string(par_value));
    else if (!strcmp(par_name, "time"))                     strcpy(shot->time, clean_string(par_value));
    else if (!strcmp(par_name, "runnumber"))                shot->runnumber = par_value_lf;
    else if (!strcmp(par_name, "matter")){
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
    else if (!strcmp(par_name, "electric_field_module"))    shot->electric_field_module = par_value_d;
    else if (!strcmp(par_name, "step_device"))              shot->step_device = par_value_d;
    else if (!strcmp(par_name, "step_host"))                shot->step_host = par_value_d;
    else if (!strcmp(par_name, "particles")){
        shot->particle_number = par_value_d;        
        shot->block_number    = par_value_d/shot->block_size+1;
    }
}

int parameter_reader(char* filename, shot_prop *shot, beam_prop *beam, run_prop *run){
    
    FILE *fp;
    char str[MAXCHAR];
    char* par_name;
    char* par_value;
 
    fp = fopen(filename, "r");
    if (fp == NULL){
        printf("Could not open file %s",filename);
        return 1;
    }
    while (fgets(str, MAXCHAR, fp) != NULL){
        if (str[0] != '#'){        
            par_name = strtok(str, "=");
            par_value = strtok(NULL, "#");
            init_taiga_props(par_name, par_value, shot, beam, run);
        }
    }    
    fclose(fp);
    strcpy(shot->name,concat(shot->shotnumber,"_",shot->time));
    return 0;
}


int runnumber_reader(char* filename, shot_prop *shot, run_prop *run){
    
    FILE *fp;
    char str[MAXCHAR];
    int runnumber;
    
    fp = fopen(filename, "r");
    if (fp == NULL){
        printf("Could not open file %s",filename);
        return 1;
    }
    fgets(str, MAXCHAR, fp);    
    sscanf(str, "%d", &runnumber);
    shot->runnumber = runnumber;
}
