__device__ double device_linear_interpolate(double *x_vector, int x_length, double *y_vector, int y_length, double x_value);
__global__ void generate_coords(taiga_globals g, taiga_commons s, beam_prop beam, beam_profile prof);
