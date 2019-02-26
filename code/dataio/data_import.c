#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_import.h"

int read_vector(double **name, char *folder, char *shotname, char *filename, bool warning_on){
    int i = 0;
    int j;
    double test;
    double *tname;
    tname = *name;
    FILE *file;
    
    char path[100];
    strcpy(path,concat(folder,"/",shotname,"/",filename));
    
    file = fopen(path,"r");
    if (file != NULL) {
        while (fscanf(file,"%lf",&test) !=EOF ) {
            i++;
        }
        fclose(file);
    
        tname = (double*)malloc(i*sizeof(double));
    
        file = fopen(path,"r");
        if (file != NULL){
            for (j=0; j<i; j++){
                fscanf(file,"%lf", &tname[j]);
            }
        }
        fclose(file);
    }else{
        if (warning_on) printf("The following file does not exists:\n%s\n\n",path);
        i = -1;
    }
        
    *name = tname;
    return i;
}

int read_vector(double **name, char *folder, char *shotname, char *filename){
    return read_vector(name, folder, shotname, filename, true);
}

int read_vector(double **name, char *folder, char *shotname, char *filename, int *successful){
    int l = read_vector(name, folder, shotname, filename, true);
    int *s;
    if (l < 0)  successful[0] = 0;
    return l;
}

int matrixColoumnReader(double **name, char *folder, char *shotname, char *filename, int coloumn_id, int total_coloumn){
    char path[100];
    strcpy(path,concat(folder,"/",shotname,"/",filename));
    return matrixColoumnReader(name, path, coloumn_id, total_coloumn);
}

int matrixColoumnReader(double **name, char *path, int coloumn_id, int total_coloumn){
    int i = 0, i2;
    int j;
    double test;
    double *tname;
    tname = *name;
    FILE *file;
    file = fopen(path,"r");
    if (file != NULL) {
        while (fscanf(file,"%lf",&test) !=EOF ) {
            i++;
        }
        
        fclose(file);
        
        i2 = i+total_coloumn-coloumn_id-1;
        tname = (double*)malloc(i2*sizeof(double));
    
        file = fopen(path,"r");
        if (file != NULL) {
            for (j=0; j<i; j++){
                if (j%total_coloumn == coloumn_id){
                    if (!fscanf(file,"%lf", &tname[j/total_coloumn]))   break;
                    printf("%lf\n",test);
                }
            }
        }        
        fclose(file);        
    }else{
        printf("The following file does not exists:\n%s\n\n",path);
        i = -1;
    }
    //*name = tname;
    return i;
}
