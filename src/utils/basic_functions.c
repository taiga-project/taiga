#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stddef.h>
#include <ctype.h>
#include "basic_functions.h"

double linear_interpolate(double *x_vector, long x_length, double *y_vector, long y_length, double x_value){
    int i;
    if (x_length != y_length)   printf("ERROR: in interpolation. Two input vectors have different length.");
    for (i=1; (i<x_length) && (x_vector[i-1]>x_value); ++i);
    if(i>1){--i;}else{i=1;}
    return y_vector[i] - (y_vector[i]-y_vector[i-1])*(x_value-x_vector[i-1])/(x_vector[i]-x_vector[i-1]);
}

char* concat(const char *s1, ...){
    const char *s;
    va_list args;
    char *r;
    size_t arg_length = 1;
    
    va_start(args, s1);
    for(s=s1; s!=NULL; s=va_arg(args, const char*)){
        arg_length += strlen(s);
    }
    va_end(args);
    r = (char*)malloc(arg_length);
    va_start(args, s1);
    for(s=s1; s!=NULL; s=va_arg(args, const char*)){
        strcat(r, s);
    }
    va_end(args);
    return r;
}

void string_to_lowercase(char* str){
    char *p;
    for(p=str; *p; ++p) *p=tolower(*p);
}