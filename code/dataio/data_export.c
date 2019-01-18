#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#if defined(_WIN32)
void mkdir(char *path, mode_t mode){
	_mkdir(path);
}
#endif

void export_data(double *dat, int Ndat, char *folder, char *runnumber, char *filename, int dat_per_line){
	char path[100];
	mkdir(folder, 0777);
	mkdir(concat(folder,"/",runnumber), 0777);
	strcpy(path,concat(folder,"/",runnumber,"/",filename));

	//! make file (open for editing)
	FILE *f = fopen(path, "a");
	if (f == NULL) {
	    printf("\nError opening file!\t\t%s\n",path);
	    exit(1);
	}
	//! write data
	for (int i=0;i<Ndat;i++){
		fprintf(f,"%le\t",dat[i]);
		if ((i+1)%dat_per_line == 0)	fprintf(f,"\n");  
	}
	fprintf(f,"\n");
	fclose(f);
}

void export_data(double *dat, int Ndat, char *folder, char *runnumber, char *filename){
	export_data(dat, Ndat, folder, runnumber, filename, INFINITY);
}

void export_data(int *dat, int Ndat, char *folder, char *runnumber, char *filename, int dat_per_line){
	char path[100];
	mkdir(folder, 0777);
	mkdir(concat(folder,"/",runnumber), 0777);
	strcpy(path,concat(folder,"/",runnumber,"/",filename));

	//! make file (open for editing)
	FILE *f = fopen(path, "a");
	if (f == NULL) {
	    printf("\nError opening file!\t\t%s\n",path);
	    exit(1);
	}

	//! write data
	for (int i=0;i<Ndat;i++){
		fprintf(f,"%d\t",dat[i]);
		if ((i+1)%dat_per_line == 0)	fprintf(f,"\n");  
	}
	fprintf(f,"\n");
	fclose(f);
}

void export_data(int *dat, int Ndat, char *folder, char *runnumber, char *filename, int dat_per_line){
	export_data(dat, Ndat, folder, runnumber, filename, INFINITY);
}

void export_header(char *dataname,char *unitname,double dat, char *folder, char *runnumber){
	char path[100];
	mkdir(folder, 0777);
	mkdir(concat(folder,"/",runnumber), 0777);
	strcpy(path,concat(folder,"/",runnumber,"/header.dat"));

	//! make file (open for editing)
	FILE *f = fopen(path, "a");
	if (f == NULL) {
	    printf("\nError opening file!\t\t%s\n",path);
	    exit(1);
	}

	//! writing to file
	fprintf(f,"%s:\t%lg %s",dataname,dat,unitname);
	fprintf(f,"\n");
	fclose(f);
}

void export_header(char *dataname, char *unitname, double dat, double dat2, char *folder, char *runnumber){
	char path[100];
	mkdir(folder, 0777);
	mkdir(concat(folder,"/",runnumber), 0777);
	strcpy(path,concat(folder,"/",runnumber,"/header.dat"));

	//! make file (open for editing)
	FILE *f = fopen(path, "a");
	if (f == NULL) {
	    printf("Error opening file!\t\t%s\n",path);
	    exit(1);
	}
	//! writing to file
	fprintf(f,"%s:\t%lg %s \t%lg %s",dataname,dat,unitname,dat2,unitname);
	fprintf(f,"\n");
	fclose(f);
}

void export_header(char *text, char *folder, char *runnumber){
	char path[100];
	mkdir(folder, 0777);
	mkdir(concat(folder,"/",runnumber), 0777);
	strcpy(path,concat(folder,"/",runnumber,"/header.dat"));

	//! make file (open for editing)
	FILE *f = fopen(path, "a");
	if (f == NULL) {
	    printf("Error opening file!\t\t%s\n",path);
	    exit(1);
	}

	//! writing to file
	fprintf(f,"%s",text);
	fprintf(f,"\n");
	fclose(f);
}
