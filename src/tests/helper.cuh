#ifndef TEST_HELPER_CUH
#define TEST_HELPER_CUH

#define LENGTH_TMP 10
void init_tmp(double *h_tmp);
void print_tmp(double *h_tmp);
void start_reference(double **h_tmp, double **d_tmp);
void end_reference(double **h_tmp, double **d_tmp);

#endif //TEST_HELPER_CUH
