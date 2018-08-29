#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*!
Load raw data
@param folder folder name
@param filename0 file name
@output double** name
*/

int vectorReader(double **name, char *folder, char *shotname, char *filename0, bool warning_on){
	int i = 0;
	int j;
	double test;
	double *tname;
	tname = *name;
	FILE *file;
	
	//! create filepath
	char filename[100];
	strcpy(filename,folder);	
	strcat(filename,"/");	
	strcat(filename,shotname);	
	strcat(filename,"/");	
	strcat(filename,filename0);	
	
	file = fopen(filename,"r");
/*		printf("Reading:\t%s\n",filename);*/
	if (file != NULL)	{
		while (fscanf(file,"%lf",&test) !=EOF )	{
			//printf("%lf\n",test);
			//printf("%d.: %lf\n\n",i,test);
			i++;
		}
		
		fclose(file);	
	
		tname = (double*)malloc(i*sizeof(double));
	
		file = fopen(filename,"r");
		if (file != NULL)	{
			for (j=0;  j<i ; j++ )	{
				fscanf(file,"%lf", &tname[j]);	
				//printf("%d %lf\n",j,tname[j]);
			}				
		}
		
		fclose(file);
		
	}else{
		if warning_on printf("The following file does not exists:\n%s\n\n",filename);
        i = -1;
	}	
	
	*name = tname;
	//printf("read %lf\n",*name[0]);
		
	return i;
	
}

int vectorReader(double **name, char *folder, char *shotname, char *filename0){
	return vectorReader(name, folder, shotname, filename0, true);
}

int vectorReader(double **name, char *folder, char *shotname, char *filename0, int *successful){
	int l = vectorReader(name, folder, shotname, filename0, true);
	int *s;
	if (l < 0){
		successful[0] = 0;
	}    
	return l;
}



/*!
Load raw data
@param folder folder name
@param filename0 file name
@output double** name
*/

int matrixColoumnReader(double **name, char *folder, char *shotname, char *filename0, int coloumn_id, int total_coloumn){
	int i = 0, i2;
	int j;
	double test;
	double *tname;
	tname = *name;
	FILE *file;
	
	//! create filepath
	char filename[100];
	strcpy(filename,folder);	
	strcat(filename,"/");	
	strcat(filename,shotname);	
	strcat(filename,"/");	
	strcat(filename,filename0);	
	
	file = fopen(filename,"r");
/*		printf("Reading:\t%s\n",filename);*/
	if (file != NULL)	{
		while (fscanf(file,"%lf",&test) !=EOF )	{
			//printf("%lf\n",test);
			//printf("%d.: %lf\n\n",i,test);
			i++;
		}
		
		fclose(file);	
		
		i2 = i+total_coloumn-coloumn_id-1;
		tname = (double*)malloc(i2*sizeof(double));
	
		file = fopen(filename,"r");
		if (file != NULL)	{
			for (j=0;  j<i ; j++ )	{
				if (j%total_coloumn == coloumn_id){
					fscanf(file,"%lf", &tname[j/total_coloumn]);	
					//printf("%d %lf\n",j,tname[j/total_coloumn]);
				}
			}				
		}
		
		fclose(file);
		
	}else{
		printf("The following file does not exists:\n%s\n\n",filename);
        i = -1;
	}	
	
	*name = tname;
	//printf("read %lf\n",*name[0]);
		
	return i;
	
}

/*!
OLD VECTORREADER
Load raw data
@param folder folder name
@param filename0 file name
@output double** name
*/


int vectorReader0(double **name, char *filename){
	int i = 0;
	int j;
	double test;
	double *tname;
	tname = *name;
	FILE *file;
	file = fopen(filename,"r");
/*		printf("Reading:\t%s\n",filename);*/
	if (file != NULL)	{
		while (fscanf(file,"%lf",&test) !=EOF )	{
			//printf("%lf\n",test);
			//printf("%d.: %lf\n\n",i,test);
			i++;
		}
		
		fclose(file);	
	
		tname = (double*)malloc(i*sizeof(double));
	
		file = fopen(filename,"r");
		if (file != NULL)	{
			for (j=0;  j<i ; j++ )	{
				fscanf(file,"%lf", &tname[j]);	
				//printf("%d %lf\n",j,tname[j]);
			}				
		}
		
		fclose(file);
		
	}else{
		printf("The following file does not exists:\n%s\n\n",filename);
        i = -1;
	}	
	
	*name = tname;
	//printf("read %lf\n",*name[0]);
		
	return i;
	
}
