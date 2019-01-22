#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXCHAR 1000

char* clean_string (char* str_in){
    char* str_out = str_in;
    str_out++;    
    str_out[strlen(str_out)-2] = 0;
    return str_out;
}

void init_taiga_props(char* par_name, char* par_value, shot_prop *shot, beam_prop *beam){
    
    double par_value_lf;    sscanf(par_value, "%lf", &par_value_lf);
    int par_value_d;        sscanf(par_value, "%d", &par_value_d);
        
    if (!strcmp(par_name, "shotnumber"))                    strcpy(shot->shotnumber, clean_string(par_value));
    else if (!strcmp(par_name, "time"))                     strcpy(shot->time, clean_string(par_value));
    else if (!strcmp(par_name, "runnumber"))                shot->runnumber = par_value_lf;
    else if (!strcmp(par_name, "matter"))                   strcpy(beam->matter, clean_string(par_value));
    else if (!strcmp(par_name, "energy"))                   beam->energy = par_value_lf;
    else if (!strcmp(par_name, "vertical_deflation"))       beam->vertical_deflation = par_value_lf;
    else if (!strcmp(par_name, "diameter"))                 beam->diameter = par_value_lf;
    else if (!strcmp(par_name, "detector"))                 strcpy(shot->detector_geometry, clean_string(par_value));
    else if (!strcmp(par_name, "detector_mask"))            strcpy(shot->detector_mask, clean_string(par_value));
    else if (!strcmp(par_name, "electric_field_module"))    shot->electric_field_module = par_value_d;
}

int parameter_reader(char* filename, shot_prop *shot, beam_prop *beam){
    
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
            init_taiga_props(par_name, par_value, shot, beam);
        }
    }    
    fclose(fp);
    strcpy(shot->name,concat(shot->shotnumber,"_",shot->time));
    beam->mass = get_mass(beam->matter);
    return 0;
}


int runnumber_reader(char* filename, shot_prop *shot){
    
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
