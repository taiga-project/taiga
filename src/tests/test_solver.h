#ifndef TEST_SOLVER_CUH
#define TEST_SOLVER_CUH

typedef struct SolverTestExtremaTag {
    double *extrema;
    int index;
    int direction;
    int counter;
}SolverTestExtrema;

double test_interp1(double x, double x0, double x1, double y0, double y1);
void init_homogeneous_field(double *X, double *B);
int get_direction(SolverTestExtrema *t, double *X, double *X_prev);
void get_extrema(SolverTestExtrema *t, double *X, double *X_prev);
double get_speed(double *X);
double run_homogeneous_field_with_solver(double timestep,
                                         void (*solve_diffeq)(double *X, double *B, double *E, double *E_prev, double eperm, double timestep));
void test_solver();

#endif //TEST_SOLVER_CUH
