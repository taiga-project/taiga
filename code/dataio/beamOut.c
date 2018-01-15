#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>


/*!
Make a data array entry to a file

@param dat array of dataset
@param Ndat number of balues in array
@param folder result folder for save
@param runnumber runnumber for data folder
@param filename0 filename

*/

void saveData1(double *dat, int Ndat, char *folder, char *runnumber, char *filename0){

	//! setting output folder
	char filename[100];
	strcpy(filename,folder);
		
	//! mkdir Windows compatible version
	#if defined(_WIN32)
		_mkdir(filename);
	#else 
		mkdir(filename, 0777); // notice that 777 is different than 0777
	#endif
		
	
	strcat(filename,"/");

//	printf("Formatted date & time : |%s|\n", buffer );

	strcat(filename,runnumber);
	strcat(filename,"/");
	
	//strcat(filename,filename0);
//	printf("fn: %s\n",filename);
	//mode_t process_mask = umask(0);
	
	
	//! mkdir Windows compatible version
	#if defined(_WIN32)
		_mkdir(filename);
	#else 
		mkdir(filename, 0777); // notice that 777 is different than 0777
	#endif
	
	
	//printf(filename);

	strcat(filename,filename0);

	//! make file
	FILE *f = fopen(filename, "w");
	if (f == NULL)
	{
	    printf("\nError opening file!\t\t%s\n",filename);
	    exit(1);
	}

	
	//! write data
	for (int i=0;i<Ndat;i++){
		fprintf(f,"%le\t",dat[i]);
	}

	fclose(f);

	
}


/*!
Add a data array entry to a file

@param dat array of dataset
@param Ndat number of balues in array
@param folder result folder for save
@param runnumber runnumber for data folder
@param filename0 filename

*/


void addData1(double *dat, int Ndat, char *folder, char *runnumber, char *filename0){

	//! setting output folder
	char filename[100];
	strcpy(filename,folder);	
	strcat(filename,"/");

	//! mkdir Windows compatible version
	#if defined(_WIN32)
		_mkdir(filename);
	#else 
		mkdir(filename, 0777); // notice that 777 is different than 0777
	#endif
//	printf("Formatted date & time : |%s|\n", buffer );


	strcat(filename,runnumber);
	strcat(filename,"/");
	
	//strcat(filename,filename0);
//	printf("fn: %s\n",filename);
	//mode_t process_mask = umask(0);
	
	//! mkdir Windows compatible version
	#if defined(_WIN32)
		_mkdir(filename);
	#else 
		mkdir(filename, 0777); // notice that 777 is different than 0777
	#endif
	
	//printf(filename);

	strcat(filename,filename0);

	//! make file (open for editing)
	FILE *f = fopen(filename, "a");
	if (f == NULL)
	{
	    printf("\nError opening file!\t\t%s\n",filename);
	    exit(1);
	}

	
	//! write data
	for (int i=0;i<Ndat;i++){
		fprintf(f,"%le\t",dat[i]);
	}
	fprintf(f,"\n");

	fclose(f);

	
}


/*!
Add a single value data entry to header.dat description file

@param dataname description for data
@param unitname unit (eg. meter) of data
@param dat data
@param folder result folder for save
@param runnumber runnumber for data folder

*/

void saveDataH(char *dataname,char *unitname,double dat, char *folder, char *runnumber){

	char filename[100];
	strcpy(filename,folder);
	
	strcat(filename,"/");
	//! mkdir Windows compatible version
	#if defined(_WIN32)
		_mkdir(filename);
	#else 
		mkdir(filename, 0777); // notice that 777 is different than 0777
	#endif
	strcat(filename,runnumber);
	strcat(filename,"/");
	
	//! mkdir Windows compatible version
	#if defined(_WIN32)
		_mkdir(filename);
	#else 
		mkdir(filename, 0777); // notice that 777 is different than 0777
	#endif
	

	strcat(filename,"header.dat");


	//! make file (open for editing)
	FILE *f = fopen(filename, "a");
	if (f == NULL)
	{
	    printf("\nError opening file!\t\t%s\n",filename);
	    exit(1);
	}

	//! writing to file
	fprintf(f,"%s:\t%lg %s",dataname,dat,unitname);
	
	
	fprintf(f,"\n");

	fclose(f);

	
}

/*!
Add a double value data entry to header.dat description file

@param dataname description for data
@param unitname unit (eg. meter) of data
@param dat first data
@param dat2 second data
@param folder result folder for save
@param runnumber runnumber for data folder

*/

void saveDataH2(char *dataname,char *unitname,double dat, double dat2, char *folder, char *runnumber){

	char filename[100];
	strcpy(filename,folder);
	
	strcat(filename,"/");

	strcat(filename,runnumber);
	strcat(filename,"/");
	
	//! mkdir Windows compatible version
	#if defined(_WIN32)
		_mkdir(filename);
	#else 
		mkdir(filename, 0777); // notice that 777 is different than 0777
	#endif
	

	strcat(filename,"header.dat");


	//! make file (open for editing)
	FILE *f = fopen(filename, "a");
	if (f == NULL)
	{
	    printf("Error opening file!\t\t%s\n",filename);
	    exit(1);
	}

	

	//! writing to file
	fprintf(f,"%s:\t%lg %s \t%lg %s",dataname,dat,unitname,dat2,unitname);
	
	fprintf(f,"\n");

	fclose(f);

	
}

/*!
Add a text entry to header.dat description file

@param text entry
@param folder result folder for save
@param runnumber runnumber for data folder

*/

void saveDataHT(char *text, char *folder, char *runnumber){

	char filename[100];
	strcpy(filename,folder);
	
	strcat(filename,"/");

	strcat(filename,runnumber);
	strcat(filename,"/");
	
	
	//! mkdir Windows compatible version
	#if defined(_WIN32)
		_mkdir(filename);
	#else 
		mkdir(filename, 0777); // notice that 777 is different than 0777
	#endif
	

	strcat(filename,"header.dat");


	//! make file (open for editing)
	FILE *f = fopen(filename, "a");
	if (f == NULL)
	{
	    printf("Error opening file!\t\t%s\n",filename);
	    exit(1);
	}

	//! writing to file

	fprintf(f,"%s",text);
	
	fprintf(f,"\n");

	fclose(f);

	
}
