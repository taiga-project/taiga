#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

double linear_interpolate(double *x_vector, long x_length, double *y_vector, long y_length, double x_value);
char* concat(const char *s1, ...);
void string_to_lowercase(char* str);

#endif //BASIC_FUNCTIONS_H
