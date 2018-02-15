// TAIGA default parameters

#define BANANA		 0		//! @param BANANA ABP: 0, banana orbits: 1
 
#define RADIONS		 1		//! @param RADIONS Real ion positions: 1, R=const 0

#define $ELM		 0		//! @param $ELM turn on <<ELM current perturbation>> mode

#define RKOLD		 0		//! @param RKOLD do not set! 0 (semi-RK: 1)

#define $3DINPUTPROF 1

#define $RENATE		0//110

#define N_BLOCKS     192		//! @param N_BLOCKS Number of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)
#define BLOCK_SIZE 	 1//30*4 		//! @param BLOCK_SIZE smaller is better (max 1M)

#define R_midions	 0.695		//! @param R_midions mid of ions at BANANA and no-RADIONS

#define $R_defl		2.3			//! radial position of deflection plates in meter -> TOROIDAL DEFLECTION

#define $deflH	 0				//! @param $deflH horizontal deflection in rad (up--down)  
#define $deflV	 0				//! @param $deflV vertical deflection in rad (left--right) -> TOROIDAL DEFLECTION


#if BANANA == 1
    #define $energy   0.5            // in keV
    #define $mass     2.013553212724 // in AMU (D)
    #define $diameter 50e-20         // in mm 
    #define $DETPOS 1 //! detector position
	#define dt		 1e-12			// timestep in seconds
	#define Nstep	 100000//00		// max step of a loop
	#define Nloop	 1000			// number of loops	
#else
	#define $energy   60				//! @param energy in keV
	#define $mass     7.016004558	//! @param mass in AMU (Li-7)
    #define $DETPOS 0.7089 //! detector position
    #define $diameter 25//4/*e-20*/      //! @param diameter in mm
    #define dt       1e-9			//! @param dt timestep in seconds
    #define Nstep    2000//000			//! @param Nstep max step of a loop
    #define Nloop    1//000				//! @param Nloop number of loops

#endif



#define ERRORCHECK() cErrorCheck(__FILE__, __LINE__)
#define PI 3.141592653589792346

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <cuda_profiler_api.h>
#include "cuda/nvToolsExt.h"

#include "dataio/fieldIn.c"

#if BANANA == 1
	#include "dataio/beamInBan.c"
#elif RADIONS == 1
	#if $3DINPUTPROF == 1
		#include "dataio/beamInFull.c"
	#elif $RENATE == 110
		#include "dataio/beamInRenate110.c"
	#else
		#include "dataio/beamIn.c"
	#endif
#else
	#include "dataio/beamInOne.c"
#endif
#include "dataio/beamOut.c"




#if RKOLD == 0
	#include "running/rk4.cu"
#else
	#include "running/rk4old.cu"
#endif

#include "running/ipol.cu"
#include "running/cyl2tor.cu"


#if RKOLD == 0
	#include "running/traj.cu"
#else
	#include "running/trajold.cu"
#endif


#include "running/ctrl.cu"

char* concat(const char *s1, const char *s2);


struct beam_prop{
    char* matter = "Li";
    double mass = 7.016004558;
    double energy = (double)$energy ;
    double diameter = (double)$diameter;
    double toroidal_deflation = (double)$deflH;   
    double vertical_deflation = (double)$deflV;
    double detector_R = (double)$DETPOS;
    
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

inline void cErrorCheck(const char *file, int line) {
  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("Error: %s\n", cudaGetErrorString(err));
    printf(" @ %s: %d\n", file, line);
    exit(-1);
  }
}

int set_cuda(){
	int num_devices, device, max_device;
	cudaGetDeviceCount(&num_devices);
	printf("Number of devices: %d\n",num_devices);
	
	if (num_devices > 1) {
        int max_multiprocessors = 0, max_device = 0;
        for (device = 0; device < num_devices; device++) {
            cudaDeviceProp properties;
            cudaGetDeviceProperties(&properties, device);
            if (max_multiprocessors < properties.multiProcessorCount) {
                max_multiprocessors = properties.multiProcessorCount;
                max_device = device;
            }
	      /*  printf("%d:%s\n",device,&properties.name);
	        printf("\tL2Cache:\t%d",	properties.l2CacheSize);
	        printf("\tNumber of cores:\t%d",	properties.warpSize);
	
	        printf("\tKernels:\t%d",	properties.concurrentKernels);
	        printf("\tThreads:\t%d",	properties.maxThreadsPerMultiProcessor);
	        printf("\tClock:\t%d",	properties.clockRate/1024);
	        printf("\n");*/
        }
        cudaSetDevice(max_device);
        for (device = 0; device < num_devices; device++) {
        	if(device==max_device) printf("-->");
            cudaDeviceProp properties;
            cudaGetDeviceProperties(&properties, device);
        	printf("\t%d:\t%s\n",device,&properties.name);
        }
        
    }
	
	cudaDeviceProp prop;
	cudaGetDevice(&max_device);
	cudaGetDeviceProperties(&prop, 0) ;  
}

double get_mass(char *s){
    double mass;
    
    if (strcmp(s,"D")==0){
        mass = 2.013553212724;
    }else if (strcmp(s,"Li")==0){
        mass = 7.016004558;
    }else if (strcmp(s,"Na")==0){
        mass = 20.0073517;
    }else if (strcmp(s,"K")==0){
        mass = 39.9639984821;
    }else if (strcmp(s,"H2")==0){
        mass = 2.013553212724;
    }else if (strcmp(s,"Li7")==0){
        mass = 7.016004558;
    }else if (strcmp(s,"Na20")==0){
        mass = 20.0073517;
    }else if (strcmp(s,"K40")==0){
        mass = 39.9639984821;
    }else{
        mass = (double)$mass;
    }
    
    return mass;

}



int spline_read_and_init(shot_prop shot, char* field_name, double ***return_s_ptr, int dimRZ){

	char* spline_folder = "input/fieldSpl";
	int suc[1] = {1};
    
	double *S0,  *s0;  vectorReader(&S0, "input/fieldSpl", shot.name, concat(field_name ,".spl11"), suc);	cudaMalloc((void **) &s0,  dimRZ); 
	double *S1,  *s1;  vectorReader(&S1, "input/fieldSpl", shot.name, concat(field_name ,".spl12"), suc);	cudaMalloc((void **) &s1,  dimRZ);
	double *S2,  *s2;  vectorReader(&S2, "input/fieldSpl", shot.name, concat(field_name ,".spl13"), suc);	cudaMalloc((void **) &s2,  dimRZ);
	double *S3,  *s3;  vectorReader(&S3, "input/fieldSpl", shot.name, concat(field_name ,".spl14"), suc);	cudaMalloc((void **) &s3,  dimRZ); 
	double *S4,  *s4;  vectorReader(&S4, "input/fieldSpl", shot.name, concat(field_name ,".spl21"), suc);	cudaMalloc((void **) &s4,  dimRZ); 
	double *S5,  *s5;  vectorReader(&S5, "input/fieldSpl", shot.name, concat(field_name ,".spl22"), suc);	cudaMalloc((void **) &s5,  dimRZ); 
	double *S6,  *s6;  vectorReader(&S6, "input/fieldSpl", shot.name, concat(field_name ,".spl23"), suc);	cudaMalloc((void **) &s6,  dimRZ); 
	double *S7,  *s7;  vectorReader(&S7, "input/fieldSpl", shot.name, concat(field_name ,".spl24"), suc);	cudaMalloc((void **) &s7,  dimRZ);
	double *S8,  *s8;  vectorReader(&S8, "input/fieldSpl", shot.name, concat(field_name ,".spl31"), suc);	cudaMalloc((void **) &s8,  dimRZ); 
	double *S9,  *s9;  vectorReader(&S9, "input/fieldSpl", shot.name, concat(field_name ,".spl32"), suc);	cudaMalloc((void **) &s9,  dimRZ); 
	double *S10, *s10; vectorReader(&S10,"input/fieldSpl", shot.name, concat(field_name ,".spl33"), suc);	cudaMalloc((void **) &s10,  dimRZ);
	double *S11, *s11; vectorReader(&S11,"input/fieldSpl", shot.name, concat(field_name ,".spl34"), suc);	cudaMalloc((void **) &s11,  dimRZ); 
	double *S12, *s12; vectorReader(&S12,"input/fieldSpl", shot.name, concat(field_name ,".spl41"), suc);	cudaMalloc((void **) &s12,  dimRZ);
	double *S13, *s13; vectorReader(&S13,"input/fieldSpl", shot.name, concat(field_name ,".spl42"), suc);	cudaMalloc((void **) &s13,  dimRZ); 
	double *S14, *s14; vectorReader(&S14,"input/fieldSpl", shot.name, concat(field_name ,".spl43"), suc);	cudaMalloc((void **) &s14,  dimRZ); 
	double *S15, *s15; vectorReader(&S15,"input/fieldSpl", shot.name, concat(field_name ,".spl44"), suc);	cudaMalloc((void **) &s15,  dimRZ);
    
	size_t dimB = 16*sizeof(double*);		
	double *S_PTR[16];	double **s_ptr;	cudaMalloc((void **) &s_ptr,  dimB);     

	S_PTR[0]  = s0; 	S_PTR[1]  = s1 ;	S_PTR[2]  = s2; 	S_PTR[3]  = s3;
	S_PTR[4]  = s4; 	S_PTR[5]  = s5; 	S_PTR[6]  = s6; 	S_PTR[7]  = s7;
	S_PTR[8]  = s8; 	S_PTR[9]  = s9; 	S_PTR[10] = s10;	S_PTR[11] = s11;
	S_PTR[12] = s12;	S_PTR[13] = s13;	S_PTR[14] = s14;	S_PTR[15] = s15;

	if (suc[0] == 1){
		cudaMemcpy(s0, S0, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s1, S1, dimRZ, cudaMemcpyHostToDevice);	
		cudaMemcpy(s2, S2, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s3, S3, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s4, S4, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s5, S5, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s6, S6, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s7, S7, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s8, S8, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s9, S9, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s10, S10, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s11, S11, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s12, S12, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s13, S13, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s14, S14, dimRZ, cudaMemcpyHostToDevice);
		cudaMemcpy(s15, S15, dimRZ, cudaMemcpyHostToDevice);			
		cudaMemcpy(s_ptr, S_PTR, dimB, cudaMemcpyHostToDevice);   
		free(S0);	free(S1);	free(S2);	free(S3);
		free(S4);	free(S5);	free(S6);	free(S7);	
		free(S8);	free(S9);	free(S10);	free(S11);	
		free(S12);	free(S13);	free(S14);	free(S15);
	}
	
	if (shot.debug == 1)
		for (int i=0;i<10;i++){
			printf("spline s0 %d %lf\n",i,S0[i]);
		}
	}
	
	*return_s_ptr = s_ptr; 
	return suc[0];	
    
}

int magnetic_field_read_and_init(shot_prop shot, double ***return_br_ptr, double ***return_bz_ptr, double ***return_bt_ptr, int dimRZ){    
    
	size_t dimB = 16*sizeof(double*);	
	double *BR_PTR[16];	double **br_ptr;	cudaMalloc((void **) &br_ptr,  dimB); 
	double *BT_PTR[16];	double **bt_ptr;	cudaMalloc((void **) &bt_ptr,  dimB); 
	double *BZ_PTR[16];	double **bz_ptr;	cudaMalloc((void **) &bz_ptr,  dimB);     

	int s;    
	s = spline_read_and_init(shot, "brad", &br_ptr, dimRZ);    
	s = spline_read_and_init(shot, "bz",   &bz_ptr, dimRZ);    
	s = spline_read_and_init(shot, "btor", &bt_ptr, dimRZ);     
    
	*return_br_ptr = br_ptr;
	*return_bz_ptr = bz_ptr;
	*return_bt_ptr = bt_ptr;

	return s;
}
int electric_field_read_and_init(shot_prop shot, double ***return_er_ptr, double ***return_ez_ptr, double ***return_et_ptr, int dimRZ){    
    
	size_t dimB = 16*sizeof(double*);	
	double *ER_PTR[16];	double **er_ptr;	cudaMalloc((void **) &er_ptr,  dimB); 
	double *ET_PTR[16];	double **et_ptr;	cudaMalloc((void **) &et_ptr,  dimB); 
	double *EZ_PTR[16];	double **ez_ptr;	cudaMalloc((void **) &ez_ptr,  dimB);     

	int s;   

	s = spline_read_and_init(shot, "erad", &er_ptr, dimRZ);    
	s = spline_read_and_init(shot, "ez",   &ez_ptr, dimRZ);    
	s = spline_read_and_init(shot, "etor", &et_ptr, dimRZ);     
    
	*return_er_ptr = er_ptr;
	*return_ez_ptr = ez_ptr;
	*return_et_ptr = et_ptr;

	return s;
}

int main(int argc, char *argv[]){
    //! @param shotname name of shot folder input folder (8714,11344,11347)	
    
    shot_prop shot;
    beam_prop beam;
	
	if (argc > 1)	shot.name = argv[1];	
	if (argc > 2)	shot.runnumber = atoi(argv[2]);
	if (argc > 3)	beam.matter = argv[3];			
	if (argc > 4)	beam.energy = atof(argv[4]);    
	if (argc > 5)	beam.vertical_deflation = atof(argv[5]);    
	if (argc > 6)	beam.diameter = atof(argv[6]);
	if (argc > 7)	beam.detector_R = atof(argv[7]);
    
	beam.mass = get_mass(beam.matter);
	printf("shotname: %s\n",shot.name);  
		
	int NX;
	int max_blocks;
	if (argc > 8)	max_blocks = atoi(argv[8])/shot.block_number+1;    
		else	max_blocks=shot.block_size;	
        
    if (argc > 9) shot.electric_field_module = atof(argv[9]);     
    
    if (argc > 10){ 
        shot.step_host = atof(argv[10]); 
        shot.step_device = 1;
    }
    
    if (argc > 11) shot.step_device = atof(argv[11]); 
	
    if (argc > 12) shot.debug = atof(argv[12]); 
    
	NX = shot.block_number * max_blocks;
	
	if ($3DINPUTPROF == 1){
        double *XR;
		NX = vectorReader0(&XR, "input/manual_profile/rad.dat");
        max_blocks = NX / shot.block_number+1;
        shot.block_size = NX;
	}	
		
	char* folder_out=concat("results/", shot.name);
	
	set_cuda();

	// set timestamp
	time_t rawtime;
	struct tm *info;
	char timestamp[80];
	sprintf(timestamp, "%d", shot.runnumber);

	// coords	
	double *X_PTR[3], **x_ptr;
	double *V_PTR[3], **v_ptr;
	size_t dimXP = 3*sizeof(double*);

	double *XR,  *xr; 
	double *XZ,  *xz;
	double *XT,  *xt;

	double *VR,  *vr; 
	double *VZ,  *vz;
	double *VT,  *vt;	
	
	printf("=============================\n");
	printf("Number of blocks (threads): %d\n", max_blocks);
	printf("Block size: %d\n", shot.block_size);
	printf("Number of particles: %d\n", NX);
	printf("Max steps on device (GPU): %d\n", shot.step_device);
	printf("Max steps on host (HDD): %d\n", shot.step_host);
	
	XR = (double*)malloc(sizeof(double)*NX);
	XZ = (double*)malloc(sizeof(double)*NX);
	XT = (double*)malloc(sizeof(double)*NX);

	VR = (double*)malloc(sizeof(double)*NX);
	VZ = (double*)malloc(sizeof(double)*NX);
	VT = (double*)malloc(sizeof(double)*NX);

	// phys. constants
	double eperm;
	eperm = 1.60217656535e-19/1.66053892173e-27/beam.mass;

	beamIn(XR, XZ, XT, VR, VZ, VT, beam.energy, eperm, NX, shot.name, beam.diameter, beam.toroidal_deflation, beam.vertical_deflation);

	//! position and velocity array allocation
	size_t dimX = NX * sizeof(double);

	cudaMalloc((void **) &xr,  dimX); 
	cudaMalloc((void **) &xz,  dimX); 
	cudaMalloc((void **) &xt,  dimX); 
	cudaMalloc((void **) &x_ptr,  dimXP); 

	cudaMalloc((void **) &vr,  dimX); 
	cudaMalloc((void **) &vz,  dimX); 
	cudaMalloc((void **) &vt,  dimX); 
	cudaMalloc((void **) &v_ptr,  dimXP); 

	//! coords pointers
	X_PTR[0] = xr;
	X_PTR[1] = xz;
	X_PTR[2] = xt;

	V_PTR[0] = vr;
	V_PTR[1] = vz;
	V_PTR[2] = vt;
	
	//! grid pointers
	double *G_PTR[2];
	double **g_ptr;
	size_t dimG = 2*sizeof(double*);	
	cudaMalloc((void **) &g_ptr,  dimG); 
	double *RG, *rg;
	double *ZG, *zg;

	// size definitions

	//! R-grid points

	int NR = vectorReader(&RG, "input/fieldSpl", shot.name, "r.spline");
	/*if ($ELM == 1){
		NR = vectorReader(&RG, "field/cuda/ipol/r.spline");
		printf("ELM mode\n");
	}else{
		NR = vectorReader(&RG, "input/fieldSpl", shot.name, "r.spline");
		printf("ELM-free mode\n");
	}*/
	size_t dimR = NR * sizeof(double);
	cudaMalloc((void **) &rg,  dimR); 
	
	//! Z-grid points
	int NZ = vectorReader(&ZG, "input/fieldSpl", shot.name, "z.spline");
	/*if ($ELM == 1){
		NZ = vectorReader(&ZG, "field/cuda/ipol/z.spline");
	}else{
		NZ = vectorReader(&ZG, "input/fieldSpl", shot.name, "z.spline");
	}*/

	size_t dimZ = NZ * sizeof(double);
	size_t dimRZ = (NR-1) * (NZ-1) * sizeof(double);
	cudaMalloc((void **) &zg,  dimZ); 

   	// grid pointer
	G_PTR[0] = rg;
	G_PTR[1] = zg;

	//! MAGN. FIELD (HOST, device) ALLOCATION          
    
    double **br_ptr, **bz_ptr, **bt_ptr;
    double **er_ptr, **ez_ptr, **et_ptr;
    
	int magnetic_field_loaded = magnetic_field_read_and_init(shot, &br_ptr,&bz_ptr,&bt_ptr, dimRZ);
	
	if (shot.electric_field_module){
		shot.electric_field_module = electric_field_read_and_init(shot, &er_ptr,&ez_ptr,&et_ptr, dimRZ);
    }
	
	// temporary test data
	double *TMP, *tmp;
	TMP = (double *)malloc(dimR);	cudaMalloc((void **) &tmp,  dimR); 

	// field direction on first
	double *BD1, *bd1;	cudaMalloc((void **) &bd1,  dimX); 	BD1 = (double *)malloc(dimX);
	double *BD2, *bd2;	cudaMalloc((void **) &bd2,  dimX);	BD2 = (double *)malloc(dimX);

	//! CUDA profiler START
	cudaProfilerStart();
	
	//! MEMCOPY (HOST2device)

	//! GRID COORDS	
	cudaMemcpy(rg, RG, dimR, cudaMemcpyHostToDevice);
	cudaMemcpy(zg, ZG, dimZ, cudaMemcpyHostToDevice);
	cudaMemcpy(g_ptr, G_PTR, dimG, cudaMemcpyHostToDevice);
	


	//! ION COORDS (HOST2device)
	cudaMemcpy(x_ptr, X_PTR, dimXP, cudaMemcpyHostToDevice);	

	//! ION SPEEDS (HOST2device)
	cudaMemcpy(v_ptr, V_PTR, dimXP, cudaMemcpyHostToDevice);
	
	// EXECUTION
	addData1(XR,NX,folder_out,timestamp,"t_rad.dat");
	addData1(XZ,NX,folder_out,timestamp,"t_z.dat");
	addData1(XZ,NX,folder_out,timestamp,"t_tor.dat");
	addData1(VR,NX,folder_out,timestamp,"t_vrad.dat");
	addData1(VZ,NX,folder_out,timestamp,"t_vz.dat");
	addData1(VZ,NX,folder_out,timestamp,"t_vtor.dat");

	//! Set CUDA timer 
	cudaEvent_t start, stop;
	float runtime;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

    if (shot.debug == 1){
        printf("ionV:  0.\t %lf\t %lf\t %lf\n",VR[0],VZ[0],VT[0]);
        printf("ionX:  0.\t %lf\t %lf\t %lf\n",XR[0],XZ[0],XT[0]);
        printf("ionX:  1.\t %lf\t %lf\t %lf\n",XR[1],XZ[1],XT[1]);
        
        printf("----------------------------------------------------------\n");
        printf("ion:  0.\t %lf\t %lf\t %lf\n",XR[0],XZ[0],XT[0]);
        printf("----------------------------------------------------------\n");
        for(int i=1; i<20; i++){
            printf("ion: %2d.\t %le\t %le\t %le\n",i,XR[i],XZ[i],XT[i]);
        }
        printf("----------------------------------------------------------\n");
    }
	// BANANA
	if (BANANA==1){
		printf("BANANA CTRL\n");
		// ION COORDS (HOST2device)
		cudaMemcpy(xr, XR, dimX, cudaMemcpyHostToDevice);
		cudaMemcpy(xz, XZ, dimX, cudaMemcpyHostToDevice);
		cudaMemcpy(xt, XT, dimX, cudaMemcpyHostToDevice);
		//cudaMemcpy(x_ptr, X_PTR, dimXP, cudaMemcpyHostToDevice);	

		// ION SPEEDS (HOST2device)
		cudaMemcpy(vr, VR, dimX, cudaMemcpyHostToDevice);
		cudaMemcpy(vz, VZ, dimX, cudaMemcpyHostToDevice);
		cudaMemcpy(vt, VT, dimX, cudaMemcpyHostToDevice);
				
		banCtrl <<< max_blocks, shot.block_size >>> (NR,NZ,br_ptr,bz_ptr,bt_ptr,g_ptr,x_ptr,bd1,bd2);
		cudaMemcpy(BD1, bd1, dimX, cudaMemcpyDeviceToHost);
		cudaMemcpy(BD2, bd2, dimX, cudaMemcpyDeviceToHost);
		addData1(BD1,NX,folder_out,timestamp,"d_b1.dat");
		addData1(BD2,NX,folder_out,timestamp,"d_b2.dat");
		addData1(VR,NX,folder_out,timestamp,"d_vr.dat");
		addData1(VZ,NX,folder_out,timestamp,"d_vz.dat");
		addData1(VT,NX,folder_out,timestamp,"d_vt.dat");
	}	
	
	for (int step_i=0;step_i<shot.step_host;step_i++){
		
		// ION COORDS (HOST2device)
		cudaMemcpy(xr, XR, dimX, cudaMemcpyHostToDevice);
		cudaMemcpy(xz, XZ, dimX, cudaMemcpyHostToDevice);
		cudaMemcpy(xt, XT, dimX, cudaMemcpyHostToDevice);
		//cudaMemcpy(x_ptr, X_PTR, dimXP, cudaMemcpyHostToDevice);	

		// ION SPEEDS (HOST2device)
		cudaMemcpy(vr, VR, dimX, cudaMemcpyHostToDevice);
		cudaMemcpy(vz, VZ, dimX, cudaMemcpyHostToDevice);
		cudaMemcpy(vt, VT, dimX, cudaMemcpyHostToDevice);
		//cudaMemcpy(v_ptr, V_PTR, dimXP, cudaMemcpyHostToDevice);
		
		
		// CUDA CODE, timer and Error catch	
		//ERRORCHECK();
		cudaEventRecord(start, 0);
		if (shot.electric_field_module){
			printf("electric_field_module ON\n");            
			ctrl <<< max_blocks, shot.block_size >>> (NR,NZ,br_ptr,bz_ptr,bt_ptr,er_ptr,ez_ptr,et_ptr,g_ptr,x_ptr,v_ptr,tmp,eperm,beam.detector_R,shot.step_device);
		}else{
			ctrl <<< max_blocks, shot.block_size >>> (NR,NZ,br_ptr,bz_ptr,bt_ptr,g_ptr,x_ptr,v_ptr,tmp,eperm,beam.detector_R,shot.step_device);
		}
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		ERRORCHECK();

		// ION COORDS (device2HOST)
		cudaMemcpy(XR, xr, dimX, cudaMemcpyDeviceToHost);
		cudaMemcpy(XZ, xz, dimX, cudaMemcpyDeviceToHost);
		cudaMemcpy(XT, xt, dimX, cudaMemcpyDeviceToHost);
		//ERRORCHECK();
		// ION SPEEDS (device2HOST)
		cudaMemcpy(VR, vr, dimX, cudaMemcpyDeviceToHost);
		cudaMemcpy(VZ, vz, dimX, cudaMemcpyDeviceToHost);
		cudaMemcpy(VT, vt, dimX, cudaMemcpyDeviceToHost);
		//ERRORCHECK();
		// Save data to files
		printf("Step\t%d/%d\n",step_i,shot.step_host);
		addData1(XR,NX,folder_out,timestamp,"t_rad.dat");
		addData1(XZ,NX,folder_out,timestamp,"t_z.dat");
		addData1(XZ,NX,folder_out,timestamp,"t_tor.dat");
		addData1(VR,NX,folder_out,timestamp,"t_vrad.dat");
		addData1(VZ,NX,folder_out,timestamp,"t_vz.dat");
		addData1(VZ,NX,folder_out,timestamp,"t_vtor.dat");
		
        if (shot.debug == 1){
            printf("Xion:  0.\t %lf\t %lf\t %lf\n",XR[0],XZ[0],XT[0]);
            printf("Xion:  1.\t %lf\t %lf\t %lf\n",XR[1],XZ[1],XT[1]);
            printf("Vion:  0.\t %lf\t %lf\t %lf\n",VR[0],VZ[0],VT[0]);
        }
        
		if (shot.step_host > 1){            
			for (int i = 1; (i < NX && XR[i] == beam.detector_R); i++){;
				if (i == NX-1) shot.step_host = step_i;
			}
        }
	}	
	
	// Get CUDA timer 
	cudaEventElapsedTime(&runtime, start, stop);
	printf ("Time for the kernel: %f s\n", runtime/1000.0);

	//! MEMCOPY (device2HOST)
	cudaMemcpy(TMP, tmp, dimR, cudaMemcpyDeviceToHost);
	if(TMP[0]!=42.24){
		printf("\n+---	-------------------+\n | Fatal error in running. | \n | The CUDA did not run well. |\n+-----------------------+\n");
	}else{
		printf("\n	Memcopy OK.\n");
	}
	
	if (shot.debug == 1)
		for (int i=0;i<10;i++) {
			printf("TMP%d\t%lf\n",i,TMP[i]);
		}
	}
	//! CUDA profiler STOP
	cudaProfilerStop();

	//! Save data to files
	saveData1(XR,NX,folder_out,timestamp,"rad.dat");
	saveData1(XZ,NX,folder_out,timestamp,"z.dat");
	saveData1(XT,NX,folder_out,timestamp,"tor.dat");
	saveData1(VR,NX,folder_out,timestamp,"vrad.dat");
	saveData1(VZ,NX,folder_out,timestamp,"vz.dat");
	saveData1(VT,NX,folder_out,timestamp,"vtor.dat");	
	
	saveDataHT(concat("Shot ID: ",shot.name),folder_out,timestamp);
	saveDataHT(concat("Run ID:  ",timestamp),folder_out,timestamp);
	saveDataHT("-----------------------------------",folder_out,timestamp);
	if(BANANA){
		saveDataHT("BANANA ORBITS",folder_out,timestamp);
	}else{		
		saveDataHT("ABP ION TRAJECTORIES",folder_out,timestamp);
		if(RADIONS){
			saveDataHT("(Real ionization position)",folder_out,timestamp); 
			if($3DINPUTPROF==1){
				saveDataHT("(3D input)",folder_out,timestamp);			
            }else if($RENATE==110){
				saveDataHT("(TS + Renate 1.1.0)",folder_out,timestamp);
			}
			
		}else{
			saveDataHT("(R=const ionization)",folder_out,timestamp);
		}
	}
	saveDataHT("-----------------------------------",folder_out,timestamp);
	saveDataH("Beam energy","keV",beam.energy,folder_out,timestamp);
	saveDataH("Atomic mass","AMU",beam.mass,folder_out,timestamp);
	saveDataH("Beam diameter","mm",beam.diameter,folder_out,timestamp);
	saveDataH2("Deflation (toroidal/vertical)","°",beam.toroidal_deflation,beam.vertical_deflation,folder_out,timestamp);
	if(!RADIONS&&!BANANA){	
		saveDataH("Ion. position (R)","m",R_midions,folder_out,timestamp);
	}
	
	saveDataH("Number of ions","",NX,folder_out,timestamp);
	saveDataHT("-----------------------------------",folder_out,timestamp);
	
	saveDataH("Detector position (R)","m",beam.detector_R,folder_out,timestamp);
	
	saveDataHT("-----------------------------------",folder_out,timestamp);
	
	saveDataH("Timestep","s",dt,folder_out,timestamp);	
	
	saveDataHT("-----------------------------------",folder_out,timestamp);
	
	saveDataH("Kernel runtime", "s", runtime/1000.0,folder_out,timestamp);
	saveDataHT("-----------------------------------",folder_out,timestamp);
	saveDataH("Number of blocks (threads)", "", max_blocks,folder_out,timestamp);
	saveDataH("Block size", "", shot.block_size,folder_out,timestamp);
	saveDataH("Length of a loop", "", shot.step_device,folder_out,timestamp);
	saveDataH("Number of loops", "", shot.step_host,folder_out,timestamp);		

	printf("\nData folder: %s/%s\n\n",folder_out,timestamp);

	//! Free CUDA
	cudaFree(x_ptr);	cudaFree(xr);	cudaFree(xz);	cudaFree(xt);
	cudaFree(g_ptr);	cudaFree(rg);	cudaFree(zg);		
	cudaFree(br_ptr);	cudaFree(bz_ptr);	cudaFree(bt_ptr);
    cudaFree(er_ptr);	cudaFree(ez_ptr);	cudaFree(et_ptr);
    /*
	cudaFree(br0);	cudaFree(br1);	cudaFree(br2);	cudaFree(br3);	
	cudaFree(br4);	cudaFree(br5);	cudaFree(br6);	cudaFree(br7);	
	cudaFree(br8);	cudaFree(br9);	cudaFree(br10);	cudaFree(br11);	
	cudaFree(br12);	cudaFree(br13);	cudaFree(br14);	cudaFree(br15);
		
	cudaFree(bz0);	cudaFree(bz1);	cudaFree(bz2);	cudaFree(bz3);
	cudaFree(bz4);	cudaFree(bz5);	cudaFree(bz6);	cudaFree(bz7);
	cudaFree(bz8);	cudaFree(bz9);	cudaFree(bz10);	cudaFree(bz11);
	cudaFree(bz12);	cudaFree(bz13);	cudaFree(bz14);	cudaFree(bz15);
	
	cudaFree(bt0);	cudaFree(bt1);	cudaFree(bt2);	cudaFree(bt3);	
	cudaFree(bt4);	cudaFree(bt5);	cudaFree(bt6);	cudaFree(bt7);	
	cudaFree(bt8);	cudaFree(bt9);	cudaFree(bt10);	cudaFree(bt11);	
	cudaFree(bt12);	cudaFree(bt13);	cudaFree(bt14);	cudaFree(bt15);	*/

	//! Free RAM
	free(RG);	free(ZG);	
	free(XR);	free(XZ);	free(XT);
	//	free(G_PTR);
	//	free(BR_PTR);	free(BZ_PTR);	free(BT_PTR);	
	
	/*free(BR0);	free(BR1);	free(BR2);	free(BR3);
	free(BR4);	free(BR5);	free(BR6);	free(BR7);	
	free(BR8);	free(BR9);	free(BR10);	free(BR11);	
	free(BR12);	free(BR13);	free(BR14);	free(BR15);
		
	free(BZ0);	free(BZ1);	free(BZ2);	free(BZ3);	
	free(BZ4);	free(BZ5);	free(BZ6);	free(BZ7);	
	free(BZ8);	free(BZ9);	free(BZ10);	free(BZ11);	
	free(BZ12);	free(BZ13);	free(BZ14);	free(BZ15);
	
	free(BT0);	free(BT1);	free(BT2);	free(BT3);
	free(BT4);	free(BT5);	free(BT6);	free(BT7);
	free(BT8);	free(BT9);	free(BT10);	free(BT11);
	free(BT12);	free(BT13);	free(BT14);	free(BT15);		*/
	
	//! FREE TMP variables (RAM, cuda)
	free(TMP);	cudaFree(tmp);

	printf("Ready.\n\n");
}

char* concat(const char *s1, const char *s2){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}
