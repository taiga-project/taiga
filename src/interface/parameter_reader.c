#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "interface/parameter_reader.h"
#include "utils/taiga_constants.h"

#define MAXCHAR 1000

char* clean_string (char* str){
    ++str;
    str[strlen(str)-2] = 0;
    return str;
}

void set_taiga_parameter(char* par_name, char* par_value, BeamProp *beam, ShotProp *shot, RunProp *run){
    double par_value_lf;    sscanf(par_value, "%lf", &par_value_lf);
    int par_value_d;        sscanf(par_value, "%d", &par_value_d);
    
    if (!strcmp(par_name, "shotnumber"))                    strcpy(shot->shotnumber, clean_string(par_value));
    else if (!strcmp(par_name, "time"))                     strcpy(shot->time, clean_string(par_value));
    else if (!strcmp(par_name, "runnumber"))                strcpy(run->runnumber, clean_string(par_value));
    else if (!strcmp(par_name, "species") || !strcmp(par_name, "beammatter") || !strcmp(par_name, "species")){
        strcpy(beam->species, clean_string(par_value));
    }
    else if (!strcmp(par_name, "charge"))                   beam->charge = par_value_lf;
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
    else if (!strcmp(par_name, "electric_field_module"))    run->is_electric_field_on = (bool)par_value_d;
    else if (!strcmp(par_name, "timestep"))                 run->timestep = par_value_lf;
    else if (!strcmp(par_name, "step_device"))              run->step_device = par_value_d;
    else if (!strcmp(par_name, "step_host"))                run->step_host = par_value_d;
    else if (!strcmp(par_name, "threads")){
        run->block_size = par_value_d;
        run->block_number    = run->particle_number/par_value_d+1;
    }
    else if (!strcmp(par_name, "particles")){
        run->particle_number = (long)par_value_d;
        run->block_number    = par_value_d/run->block_size+1;
    }
    else if (!strcmp(par_name, "solver"))                   set_solver(run, clean_string(par_value));
    else if (!strcmp(par_name, "field_interpolation_method")){
        set_field_interpolation_method(run, clean_string(par_value));
    }
    else if (!strcmp(par_name, "secondary_ionisation"))     run->is_ionisation_on = (bool)par_value_d;
    else if (!strcmp(par_name, "electric_field_value") || !strcmp(par_name, "matlab"))
        ;
    else{
        printf("Warning: Undefined parameter:\n\t%s: %s\n", par_name, par_value);
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
            set_taiga_parameter(par_name, par_value, beam, shot, run);
        }
    }
    fclose(fp);
    strcpy(shot->name, concat(shot->shotnumber, "_", shot->time, NULL));

    strcpy(shot->long_name, concat(shot->shotnumber, "_", shot->time, "/", beam->species, "/", to_str(beam->energy), NULL));
    strcpy(run->folder_out, concat("results/", shot->shotnumber, "_", shot->time, NULL));
    return 0;
}

int runnumber_reader(ShotProp *shot, RunProp *run){
    if (!strcmp(run->runnumber, UNDEFINED_RUNNUMBER)){
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
    return 0;
}

void set_solver(RunProp *run, char* solver){
    string_to_lowercase(solver);
    if (!strcmp(solver, "rk") || !strcmp(solver, "rk45") || !strcmp(solver, "runge-kutta") || !strcmp(solver, "runge--kutta")) {
        run->solver = SOLVER_RK45;
    }else if (!strcmp(solver, "rkn") || !strcmp(solver, "nystrom") || !strcmp(solver, "runge-kutta-nystrom") || !strcmp(solver, "runge--kutta--nystrom")){
        run->solver = SOLVER_RUNGE_KUTTA_NYSTROM;
    }else if (!strcmp(solver, "verlet")){
        run->solver = SOLVER_VERLET;
    }else if (!strcmp(solver, "yoshida")){
        run->solver = SOLVER_YOSHIDA;
    }else{
        printf("Warning: Unvalid numerical solver: %s\nSet to Runge--Kutta solver [rk45] .", solver);
    }
}

void set_field_interpolation_method(RunProp *run, char* method){
    string_to_lowercase(method);
    if (!strcmp(method, "spline")){
        run->field_interpolation_method = CUBIC_SPLINE;
    }else if (!strcmp(method, "bspline") || !strcmp(method, "b-spline")){
        run->field_interpolation_method = CUBIC_BSPLINE;
    }else{
        printf("Warning: Unvalid interpolation method: %s\nSet to bicubic spline [spline].", method);
    }
}
