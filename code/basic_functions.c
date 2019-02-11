char* concat(const char *s1, const char *s2){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

char* concat(const char *s1, const char *s2, const char *s3){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+1);
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    return result;
}

char* concat(const char *s1, const char *s2, const char *s3, const char *s4){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+strlen(s4)+1);
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    strcat(result, s4);
    return result;
}

char* concat(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+strlen(s3)+strlen(s4)+strlen(s5)+1);
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    strcat(result, s4);
    strcat(result, s5);
    return result;
}

double linear_interpolate(double *x_vector, int x_length, double *y_vector, int y_length, double x_value){
    int i;
    if (x_length != y_length)   printf("ERROR: in interpolation. Two input vectors have different length.");    
    for (i=1; (i<x_length) && (x_vector[i-1]>x_value); i++);    
    if(i>1){--i;}else{i=1;}    
    return y_vector[i] - (y_vector[i]-y_vector[i-1])*(x_value-x_vector[i-1])/(x_vector[i]-x_vector[i-1]);
}

int get_array_size(double *array){
    return (int)(sizeof(array)/sizeof(array[0]));
} 
