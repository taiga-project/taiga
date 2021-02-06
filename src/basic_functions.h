void set_cuda(int debug_flag);
inline void cErrorCheck(const char *file, int line);

int get_array_size(double *array);
double linear_interpolate(double *x_vector, long x_length, double *y_vector, long y_length, double x_value);

void init_dir(char *folder, char *runnumber);
void init_dir(char *folder, char *runnumber, char *subdir);

char* concat(const char *s1, ...);

