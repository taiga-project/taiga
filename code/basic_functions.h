void set_cuda(int debug_flag);
inline void cErrorCheck(const char *file, int line);

int get_array_size(double *array);
double linear_interpolate(double *x_vector, int x_length, double *y_vector, int y_length, double x_value);

void init_dir(char *folder, char *runnumber);
void init_dir(char *folder, char *runnumber, char *subdir);

char* concat(const char *s1, const char *s2);
char* concat(const char *s1, const char *s2, const char *s3);
char* concat(const char *s1, const char *s2, const char *s3, const char *s4);
char* concat(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5);
char* concat(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5, const char *s6);
char* concat(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5, const char *s6, const char *s7);
