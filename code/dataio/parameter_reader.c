#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXCHAR 1000

void init_taiga_props(char* par_name, char* par_value, shot_prop *shot, beam_prop *beam){
    
    double par_value_lf;
    char par_value_s[100];
    sscanf(par_value, "%lf", &par_value_lf);
        
    if (!strcmp(par_name, "shotnumber"))                {
        strcpy(shot->shotnumber, par_value);
        printf("-->%s\n",shot->shotnumber);        
    }
    else if (!strcmp(par_name, "time"))                 strcpy(shot->time, par_value);
    else if (!strcmp(par_name, "runnumber"))            shot->runnumber = par_value_lf;
    else if (!strcmp(par_name, "matter"))               strcpy(beam->matter, par_value);
    else if (!strcmp(par_name, "energy"))               beam->energy = par_value_lf;
    else if (!strcmp(par_name, "vertical_deflation"))   beam->vertical_deflation = par_value_lf;
    else if (!strcmp(par_name, "diameter"))             beam->diameter = par_value_lf;    
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
            printf("%s\t%s\n",par_name,par_value);
            init_taiga_props(par_name, par_value, shot, beam);
            printf("==>%s\n",shot->shotnumber);
        }
    }    
    fclose(fp);
    printf("DATA:%s\t%s\t%lf\n",shot->shotnumber, shot->time, beam->energy);
    printf("shot->name\t%s\n",concat(shot->shotnumber,"_",shot->time));
    shot->name = concat(shot->shotnumber,"_",shot->time);
    beam->mass = get_mass(beam->matter);
    return 0;
}
