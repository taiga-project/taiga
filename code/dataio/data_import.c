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

int read_matrix_column(double **name, char *path, int coloumn_id){
    int i, j;

    char tmp[1000];
    char *token;

    double value;
    double *tname;
    tname = *name;
    FILE *file;

    file = fopen(path,"r");
    if (file != NULL) {
        for (i=0; fgets(tmp, sizeof tmp, file) != NULL; i++)   ;
        fclose(file);
        tname = (double*)malloc(i*sizeof(double));

        file = fopen(path,"r");
        for (i=0; (fgets(tmp, sizeof tmp, file) != NULL); NULL){
            token = strtok( tmp, " " );
            for(j = 0; j<coloumn_id && token != NULL; j++){
                if (j == coloumn_id-1){
                    value = atof(token);
                    if (value || (token[0]=='0') || (strncmp(token,"-0",2) == 0) ){
                        tname[i] = value;
                        i++;
                    }
                }
                token = strtok( NULL, " " );
            }
            if (isnan(value)){
                printf("NAN DATA:n %s\n", tmp);
            }
        }

        fclose(file);
    }else{
        printf("The following file does not exists:\n%s\n\n",path);
        i = -1;
    }
    *name = tname;
    return i;
}
