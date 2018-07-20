struct beam_prop{
	char* matter = "Li";
	double mass = 7.016004558;
	double energy = (double)$energy ;
	double diameter = (double)$diameter;
	double toroidal_deflation = (double)$deflH;   
	double vertical_deflation = (double)$deflV;
	
};

struct shot_prop{
	char* name = "11347";
	int runnumber = 0;  
	int electric_field_module = 0;
	int debug = 0;
	int block_size = BLOCK_SIZE;
	int block_number = N_BLOCKS;
	int step_host = 1; // on HDD
	int step_device = 2000; // on GPU
};

int set_cuda();
inline void cErrorCheck(const char *file, int line);

int main(int argc, char *argv[]);

char* concat(const char *s1, const char *s2);
double get_mass(char *s);
void fill_detector(double *DETECTOR, char* values);

int spline_read_and_init(shot_prop shot, char* field_name, double ***return_s_ptr, int dimRZ);

int magnetic_field_read_and_init(shot_prop shot, double ***return_br_ptr, double ***return_bz_ptr, double ***return_bt_ptr, int dimRZ);
int electric_field_read_and_init(shot_prop shot, double ***return_er_ptr, double ***return_ez_ptr, double ***return_et_ptr, int dimRZ);

void debug_message_init(double* XR, double* XZ, double* XT, double* VR, double* VZ, double* VT);
void debug_message_run(double* XR, double* XZ, double* XT, double* VR, double* VZ, double* VT);
