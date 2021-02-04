__device__ double device_linear_interpolate(double *x_vector, double *y_vector, long length, double x_value);
__global__ void generate_coords(TaigaGlobals *globals, BeamProp beam, BeamProfile *prof);

