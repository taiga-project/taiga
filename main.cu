// TAIGA default parameters

#define BANANA		 0		//! @param BANANA ABP: 0, banana orbits: 1
 
#define RADIONS		 1		//! @param RADIONS Real ion positions: 1, R=const 0

#define $ELM		 0		//! @param $ELM turn on <<ELM current perturbation>> mode

#define RKOLD		 0		//! @param RKOLD do not set! 0 (semi-RK: 1)

#define $FASTMODE    0//1  	//! @param $FASTMODE fastmode do not set! 0

#define $RENATE		110

#define N_BLOCKS     192		//! @param N_BLOCKS Number of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)
#define BLOCK_SIZE 	 1//30*4 		//! @param BLOCK_SIZE smaller is better (max 1M)

#define R_midions	 0.695		//! @param R_midions mid of ions at BANANA and no-RADIONS

#define $R_defl		2.3			//! radial position of deflection plates in meter -> TOROIDAL DEFLECTION

#define $deflH	 0				//! @param $deflH horizontal deflection in rad (up--down)  
#define $deflV	 0				//! @param $deflV vertical deflection in rad (left--right) -> TOROIDAL DEFLECTION

#define $DETPOS 0.7089 //! detector position

#if BANANA == 1
    #define $energy   0.5            // in keV
    #define $mass     2.013553212724 // in AMU (D)
    #define $diameter 50e-20         // in mm 
	#define dt		 1e-12			// timestep in seconds
	#define Nstep	 100000//00		// max step of a loop
	#define Nloop	 1000			// number of loops	
#else
	#define $energy   60				//! @param energy in keV
	#define $mass     7.016004558	//! @param mass in AMU (Li-7)
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
	#if $RENATE == 110
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


inline void cErrorCheck(const char *file, int line) {
  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("Error: %s\n", cudaGetErrorString(err));
    printf(" @ %s: %d\n", file, line);
    exit(-1);
  }
}

double get_mass(char *s){
    double mass;
    
    if (strcmp(s,"D")){
        mass = 2.013553212724;
    }else if (strcmp(s,"Li")){
        mass = 7.016004558;
    }else if (strcmp(s,"Na")){
        mass = 20.0073517;
    }else if (strcmp(s,"K")){
        mass = 39.9639984821;
    }else if (strcmp(s,"H2")){
        mass = 2.013553212724;
    }else if (strcmp(s,"Li7")){
        mass = 7.016004558;
    }else if (strcmp(s,"Na20")){
        mass = 20.0073517;
    }else if (strcmp(s,"K40")){
        mass = 39.9639984821;
    }else{
        mass = (double)$mass;
    }
    
    return mass;

}

int main(int argc, char *argv[]){
    //! @param shotname name of shot folder input folder (8714,11344,11347)	
	char* shotname;
	char *beammatter;
	if (argc >= 2){
		shotname = argv[1];	
	}else{
		shotname = "11347";
	}	
	
	printf("shotname: %s\n",shotname);
	
	
	if (argc >= 3){
		beammatter = argv[2];	
	}else{
		beammatter = "Li";
	}	
	double mass = get_mass(beammatter);
	
	
	double energy=(double)$energy;
	if (argc >= 4){
		energy = atof(argv[3]);
    }
	
	double deflV=(double)$deflV;
    double deflH=(double)$deflH;
	if (argc >= 5){
		deflV = atof(argv[4]);
    }
    
    double diameter=(double)$diameter;
	if (argc >= 6){
		diameter = atof(argv[5]);
    }
		
	int NX;
	int max_blocks;
	if (argc >= 7){
		max_blocks = atoi(argv[6])/BLOCK_SIZE+1;
		//printf("max blocks: %d\n\n",max_blocks);
        //NX = atoi(argv[1]); //for the future
    }else{        
        //NX = BLOCK_SIZE*N_BLOCKS;
        max_blocks=N_BLOCKS;
	}
	
    NX = BLOCK_SIZE * max_blocks;
    

    
		
	char* folder_out=concat("results/", shotname);//! io properties folder
	// card settings
	
	
	
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
	cudaGetDeviceProperties( &prop, 0) ;
 /*0*/	
//	int BLOCK_SIZE = 1;//prop.maxThreadsPerBlock;
//	if(BLOCK_SIZE<1) BLOCK_SIZE=1;
//	int N_BLOCKS = 192;//0;



	// detector position
	double l_ri;

	// FOR ABP
	if (!BANANA) {
		l_ri = $DETPOS;
	// FOR BANANA ORBITS
	}else{
		l_ri = 1.0000; // just for taiga does not stop, there is no physical meaning
	}


	// phys. constants
	double eperm;

	// set timestamp

	time_t rawtime;
	struct tm *info;
	char timestamp[80];
	time( &rawtime );
  	info = localtime( &rawtime );
  	strftime(timestamp,80,"%d%b%Y_%H%M%S", info);
	
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
	printf("Block size: %d\n", BLOCK_SIZE);
	printf("Number of parts: %d\n", NX);
	printf("length of a loop: %d\n", Nstep);
	printf("Number of loops: %d\n", Nloop);


	
	XR = (double*)malloc(sizeof(double)*NX);
	XZ = (double*)malloc(sizeof(double)*NX);
	XT = (double*)malloc(sizeof(double)*NX);

	VR = (double*)malloc(sizeof(double)*NX);
	VZ = (double*)malloc(sizeof(double)*NX);
	VT = (double*)malloc(sizeof(double)*NX);


	//time_t t;

	eperm = 1.60217656535e-19/1.66053892173e-27/mass;

	beamIn(XR, XZ, XT, VR, VZ, VT, energy, eperm, NX, shotname, diameter, deflH, deflV);
	/*XR[0] = 0.72;
	XZ[0] = 0.00;
	XT[0] = 0.00;*/
	

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

	int NR = vectorReader(&RG, "input/fieldSpl", shotname, "r.spline");
	/*if ($ELM == 1){
		NR = vectorReader(&RG, "field/cuda/ipol/r.spline");
		printf("ELM mode\n");
	}else{
		NR = vectorReader(&RG, "input/fieldSpl", shotname, "r.spline");
		printf("ELM-free mode\n");
	}*/
	size_t dimR = NR * sizeof(double);
	cudaMalloc((void **) &rg,  dimR); 
	
	//! Z-grid points
	int NZ = vectorReader(&ZG, "input/fieldSpl", shotname, "z.spline");
	/*if ($ELM == 1){
		NZ = vectorReader(&ZG, "field/cuda/ipol/z.spline");
	}else{
		NZ = vectorReader(&ZG, "input/fieldSpl", shotname, "z.spline");
	}*/
	//int NRZ = NR*NZ;
	size_t dimZ = NZ * sizeof(double);
	size_t dimRZ = (NR-1) * (NZ-1) * sizeof(double);
	cudaMalloc((void **) &zg,  dimZ); 

   	// grid pointer
	G_PTR[0] = rg;
	G_PTR[1] = zg;


	//! MAGN. FIELD (HOST, device) ALLOCATION
	
	//!rad
	double *BR0,  *br0;  vectorReader(&BR0, "input/fieldSpl", shotname, "brad.spl11");	cudaMalloc((void **) &br0,  dimRZ); 
	double *BR1,  *br1;  vectorReader(&BR1, "input/fieldSpl", shotname, "brad.spl12");	cudaMalloc((void **) &br1,  dimRZ);
	double *BR2,  *br2;  vectorReader(&BR2, "input/fieldSpl", shotname, "brad.spl13");	cudaMalloc((void **) &br2,  dimRZ);
	double *BR3,  *br3;  vectorReader(&BR3, "input/fieldSpl", shotname, "brad.spl14");	cudaMalloc((void **) &br3,  dimRZ); 
	double *BR4,  *br4;  vectorReader(&BR4, "input/fieldSpl", shotname, "brad.spl21");	cudaMalloc((void **) &br4,  dimRZ); 
	double *BR5,  *br5;  vectorReader(&BR5, "input/fieldSpl", shotname, "brad.spl22");	cudaMalloc((void **) &br5,  dimRZ); 
	double *BR6,  *br6;  vectorReader(&BR6, "input/fieldSpl", shotname, "brad.spl23");	cudaMalloc((void **) &br6,  dimRZ); 
	double *BR7,  *br7;  vectorReader(&BR7, "input/fieldSpl", shotname, "brad.spl24");	cudaMalloc((void **) &br7,  dimRZ);
	double *BR8,  *br8;  vectorReader(&BR8, "input/fieldSpl", shotname, "brad.spl31");	cudaMalloc((void **) &br8,  dimRZ); 
	double *BR9,  *br9;  vectorReader(&BR9, "input/fieldSpl", shotname, "brad.spl32");	cudaMalloc((void **) &br9,  dimRZ); 
	double *BR10, *br10; vectorReader(&BR10,"input/fieldSpl", shotname, "brad.spl33");	cudaMalloc((void **) &br10,  dimRZ);
	double *BR11, *br11; vectorReader(&BR11,"input/fieldSpl", shotname, "brad.spl34");	cudaMalloc((void **) &br11,  dimRZ); 
	double *BR12, *br12; vectorReader(&BR12,"input/fieldSpl", shotname, "brad.spl41");	cudaMalloc((void **) &br12,  dimRZ);
	double *BR13, *br13; vectorReader(&BR13,"input/fieldSpl", shotname, "brad.spl42");	cudaMalloc((void **) &br13,  dimRZ); 
	double *BR14, *br14; vectorReader(&BR14,"input/fieldSpl", shotname, "brad.spl43");	cudaMalloc((void **) &br14,  dimRZ); 
	double *BR15, *br15; vectorReader(&BR15,"input/fieldSpl", shotname, "brad.spl44");	cudaMalloc((void **) &br15,  dimRZ);

	//!tor
	double *BT0,  *bt0;  vectorReader(&BT0, "input/fieldSpl", shotname, "btor.spl11");	cudaMalloc((void **) &bt0,  dimRZ); 
	double *BT1,  *bt1;  vectorReader(&BT1, "input/fieldSpl", shotname, "btor.spl12");	cudaMalloc((void **) &bt1,  dimRZ); 
	double *BT2,  *bt2;  vectorReader(&BT2, "input/fieldSpl", shotname, "btor.spl13");	cudaMalloc((void **) &bt2,  dimRZ); 
	double *BT3,  *bt3;  vectorReader(&BT3, "input/fieldSpl", shotname, "btor.spl14");	cudaMalloc((void **) &bt3,  dimRZ); 
	double *BT4,  *bt4;  vectorReader(&BT4, "input/fieldSpl", shotname, "btor.spl21");	cudaMalloc((void **) &bt4,  dimRZ); 
	double *BT5,  *bt5;  vectorReader(&BT5, "input/fieldSpl", shotname, "btor.spl22");	cudaMalloc((void **) &bt5,  dimRZ); 
	double *BT6,  *bt6;  vectorReader(&BT6, "input/fieldSpl", shotname, "btor.spl23");	cudaMalloc((void **) &bt6,  dimRZ); 
	double *BT7,  *bt7;  vectorReader(&BT7, "input/fieldSpl", shotname, "btor.spl24");	cudaMalloc((void **) &bt7,  dimRZ);
	double *BT8,  *bt8;  vectorReader(&BT8, "input/fieldSpl", shotname, "btor.spl31");	cudaMalloc((void **) &bt8,  dimRZ); 
	double *BT9,  *bt9;  vectorReader(&BT9, "input/fieldSpl", shotname, "btor.spl32");	cudaMalloc((void **) &bt9,  dimRZ); 
	double *BT10, *bt10; vectorReader(&BT10,"input/fieldSpl", shotname, "btor.spl33");	cudaMalloc((void **) &bt10,  dimRZ); 
	double *BT11, *bt11; vectorReader(&BT11,"input/fieldSpl", shotname, "btor.spl34");	cudaMalloc((void **) &bt11,  dimRZ); 
	double *BT12, *bt12; vectorReader(&BT12,"input/fieldSpl", shotname, "btor.spl41");	cudaMalloc((void **) &bt12,  dimRZ); 
	double *BT13, *bt13; vectorReader(&BT13,"input/fieldSpl", shotname, "btor.spl42");	cudaMalloc((void **) &bt13,  dimRZ);
	double *BT14, *bt14; vectorReader(&BT14,"input/fieldSpl", shotname, "btor.spl43");	cudaMalloc((void **) &bt14,  dimRZ); 
	double *BT15, *bt15; vectorReader(&BT15,"input/fieldSpl", shotname, "btor.spl44");	cudaMalloc((void **) &bt15,  dimRZ);
	
	//!z
	double *BZ0,  *bz0;  vectorReader(&BZ0, "input/fieldSpl", shotname, "bz.spl11");	cudaMalloc((void **) &bz0,  dimRZ); 
	double *BZ1,  *bz1;  vectorReader(&BZ1, "input/fieldSpl", shotname, "bz.spl12");	cudaMalloc((void **) &bz1,  dimRZ); 
	double *BZ2,  *bz2;  vectorReader(&BZ2, "input/fieldSpl", shotname, "bz.spl13");	cudaMalloc((void **) &bz2,  dimRZ); 
	double *BZ3,  *bz3;  vectorReader(&BZ3, "input/fieldSpl", shotname, "bz.spl14");	cudaMalloc((void **) &bz3,  dimRZ); 
	double *BZ4,  *bz4;  vectorReader(&BZ4, "input/fieldSpl", shotname, "bz.spl21");	cudaMalloc((void **) &bz4,  dimRZ);
	double *BZ5,  *bz5;  vectorReader(&BZ5, "input/fieldSpl", shotname, "bz.spl22");	cudaMalloc((void **) &bz5,  dimRZ); 
	double *BZ6,  *bz6;  vectorReader(&BZ6, "input/fieldSpl", shotname, "bz.spl23");	cudaMalloc((void **) &bz6,  dimRZ);
	double *BZ7,  *bz7;  vectorReader(&BZ7, "input/fieldSpl", shotname, "bz.spl24");	cudaMalloc((void **) &bz7,  dimRZ);
	double *BZ8,  *bz8;  vectorReader(&BZ8, "input/fieldSpl", shotname, "bz.spl31");	cudaMalloc((void **) &bz8,  dimRZ); 
	double *BZ9,  *bz9;  vectorReader(&BZ9, "input/fieldSpl", shotname, "bz.spl32");	cudaMalloc((void **) &bz9,  dimRZ); 
	double *BZ10, *bz10; vectorReader(&BZ10,"input/fieldSpl", shotname, "bz.spl33");	cudaMalloc((void **) &bz10,  dimRZ); 
	double *BZ11, *bz11; vectorReader(&BZ11,"input/fieldSpl", shotname, "bz.spl34");	cudaMalloc((void **) &bz11,  dimRZ); 
	double *BZ12, *bz12; vectorReader(&BZ12,"input/fieldSpl", shotname, "bz.spl41");	cudaMalloc((void **) &bz12,  dimRZ);
	double *BZ13, *bz13; vectorReader(&BZ13,"input/fieldSpl", shotname, "bz.spl42");	cudaMalloc((void **) &bz13,  dimRZ); 
	double *BZ14, *bz14; vectorReader(&BZ14,"input/fieldSpl", shotname, "bz.spl43");	cudaMalloc((void **) &bz14,  dimRZ);
	double *BZ15, *bz15; vectorReader(&BZ15,"input/fieldSpl", shotname, "bz.spl44");	cudaMalloc((void **) &bz15,  dimRZ);
	

	// magnetic field pointer array
	// magnetic field (HOST, device)
	size_t dimB = 16*sizeof(double*);
		
	double *BR_PTR[16];	double **br_ptr;	
	cudaMalloc((void **) &br_ptr,  dimB); 
	double *BT_PTR[16];	double **bt_ptr;
	cudaMalloc((void **) &bt_ptr,  dimB); 
	double *BZ_PTR[16];	double **bz_ptr;
	cudaMalloc((void **) &bz_ptr,  dimB); 	
		
	
	
	//! MAGN. FIELD POINTERS
	
	//!rad	
	BR_PTR[0] = br0;	BR_PTR[1] = br1;	BR_PTR[2] = br2;	BR_PTR[3] = br3;
	BR_PTR[4] = br4;	BR_PTR[5] = br5;	BR_PTR[6] = br6;	BR_PTR[7] = br7;
	BR_PTR[8] = br8;	BR_PTR[9] = br9;	BR_PTR[10] = br10;	BR_PTR[11] = br11;
	BR_PTR[12] = br12;	BR_PTR[13] = br13;	BR_PTR[14] = br14;	BR_PTR[15] = br15;
	
	//!tor
	BT_PTR[0] = bt0;	BT_PTR[1] = bt1;	BT_PTR[2] = bt2;	BT_PTR[3] = bt3;
	BT_PTR[4] = bt4;	BT_PTR[5] = bt5;	BT_PTR[6] = bt6;	BT_PTR[7] = bt7;
	BT_PTR[8] = bt8;	BT_PTR[9] = bt9;	BT_PTR[10] = bt10;	BT_PTR[11] = bt11;
	BT_PTR[12] = bt12;	BT_PTR[13] = bt13;	BT_PTR[14] = bt14;	BT_PTR[15] = bt15;
	
	//!z
	BZ_PTR[0] = bz0;	BZ_PTR[1] = bz1;	BZ_PTR[2] = bz2;	BZ_PTR[3] = bz3;
	BZ_PTR[4] = bz4;	BZ_PTR[5] = bz5;	BZ_PTR[6] = bz6;	BZ_PTR[7] = bz7;
	BZ_PTR[8] = bz8;	BZ_PTR[9] = bz9;	BZ_PTR[10] = bz10;	BZ_PTR[11] = bz11;
	BZ_PTR[12] = bz12;	BZ_PTR[13] = bz13;	BZ_PTR[14] = bz14;	BZ_PTR[15] = bz15;
	
	
	
		
	// temporary test data
	double *TMP, *tmp;
	TMP = (double *)malloc(dimR);
	cudaMalloc((void **) &tmp,  dimR); 
	
	
		
	// field direction on first
	double *BD1, *bd1;
	double *BD2, *bd2;
	BD1 = (double *)malloc(dimX);
	BD2 = (double *)malloc(dimX);
	cudaMalloc((void **) &bd1,  dimX); 
	cudaMalloc((void **) &bd2,  dimX); 
	
	

	//! CUDA profiler START
	cudaProfilerStart();
	
	//! MEMCOPY (HOST2device)
	
	


	//! GRID COORDS	
	cudaMemcpy(rg, RG, dimR, cudaMemcpyHostToDevice);
	cudaMemcpy(zg, ZG, dimZ, cudaMemcpyHostToDevice);
	cudaMemcpy(g_ptr, G_PTR, dimG, cudaMemcpyHostToDevice);
	
	//! MAGNETIC FIELD
	
	//!rad	
	cudaMemcpy(br0, BR0, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br1, BR1, dimRZ, cudaMemcpyHostToDevice);	
	cudaMemcpy(br2, BR2, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br3, BR3, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br4, BR4, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br5, BR5, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br6, BR6, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br7, BR7, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br8, BR8, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br9, BR9, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br10, BR10, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br11, BR11, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br12, BR12, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br13, BR13, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br14, BR14, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(br15, BR15, dimRZ, cudaMemcpyHostToDevice);			
	cudaMemcpy(br_ptr, BR_PTR, dimB, cudaMemcpyHostToDevice);
	
	//!tor
	cudaMemcpy(bt0, BT0, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt1, BT1, dimRZ, cudaMemcpyHostToDevice);	
	cudaMemcpy(bt2, BT2, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt3, BT3, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt4, BT4, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt5, BT5, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt6, BT6, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt7, BT7, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt8, BT8, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt9, BT9, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt10, BT10, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt11, BT11, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt12, BT12, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt13, BT13, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt14, BT14, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bt15, BT15, dimRZ, cudaMemcpyHostToDevice);			
	cudaMemcpy(bt_ptr, BT_PTR, dimB, cudaMemcpyHostToDevice);
	
	//!z
	cudaMemcpy(bz0, BZ0, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz1, BZ1, dimRZ, cudaMemcpyHostToDevice);	
	cudaMemcpy(bz2, BZ2, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz3, BZ3, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz4, BZ4, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz5, BZ5, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz6, BZ6, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz7, BZ7, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz8, BZ8, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz9, BZ9, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz10, BZ10, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz11, BZ11, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz12, BZ12, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz13, BZ13, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz14, BZ14, dimRZ, cudaMemcpyHostToDevice);
	cudaMemcpy(bz15, BZ15, dimRZ, cudaMemcpyHostToDevice);			
	cudaMemcpy(bz_ptr, BZ_PTR, dimB, cudaMemcpyHostToDevice);



	//! ION COORDS (HOST2device)
	/*cudaMemcpy(xr, XR, dimX, cudaMemcpyHostToDevice);
	cudaMemcpy(xz, XZ, dimX, cudaMemcpyHostToDevice);
	cudaMemcpy(xt, XT, dimX, cudaMemcpyHostToDevice);*/
	cudaMemcpy(x_ptr, X_PTR, dimXP, cudaMemcpyHostToDevice);	

	//! ION SPEEDS (HOST2device)
	/*cudaMemcpy(vr, VR, dimX, cudaMemcpyHostToDevice);
	cudaMemcpy(vz, VZ, dimX, cudaMemcpyHostToDevice);
	cudaMemcpy(vt, VT, dimX, cudaMemcpyHostToDevice);*/
	cudaMemcpy(v_ptr, V_PTR, dimXP, cudaMemcpyHostToDevice);

	
	// EXECUTION	




	addData1(XR,NX,folder_out,timestamp,"t_rad.dat");
	addData1(XZ,NX,folder_out,timestamp,"t_z.dat");
	addData1(XZ,NX,folder_out,timestamp,"t_tor.dat");

	//! Set CUDA timer 
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	printf("ionV:  0.\t %lf\t %lf\t %lf\n",VR[0],VZ[0],VT[0]);
	printf("ionX:  0.\t %lf\t %lf\t %lf\n",XR[0],XZ[0],XT[0]);
	printf("ionX:  1.\t %lf\t %lf\t %lf\n",XR[1],XZ[1],XT[1]);
	
	
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
				
		banCtrl <<< max_blocks, BLOCK_SIZE >>> (NR,NZ,br_ptr,bz_ptr,bt_ptr,g_ptr,x_ptr,bd1,bd2);
		cudaMemcpy(BD1, bd1, dimX, cudaMemcpyDeviceToHost);
		cudaMemcpy(BD2, bd2, dimX, cudaMemcpyDeviceToHost);
		addData1(BD1,NX,folder_out,timestamp,"d_b1.dat");
		addData1(BD2,NX,folder_out,timestamp,"d_b2.dat");
		addData1(VR,NX,folder_out,timestamp,"d_vr.dat");
		addData1(VZ,NX,folder_out,timestamp,"d_vz.dat");
		addData1(VT,NX,folder_out,timestamp,"d_vt.dat");
	}	
	
	for (int step_i=0;step_i<Nloop;step_i++){

		
		if ($FASTMODE==0){
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
		}else{
			double *PROFR, *PROFD;
			int N_beamdens = vectorReader0(&PROFR,"dataio/data/prof_r.dat");
			vectorReader0(&PROFD,"dataio/data/prof_d.dat");
		}
		
		// CUDA CODE, timer and Error catch	
		//ERRORCHECK();
		cudaEventRecord(start, 0);
		ctrl <<< max_blocks, BLOCK_SIZE >>> (NR,NZ,br_ptr,bz_ptr,bt_ptr,g_ptr,x_ptr,v_ptr,tmp,eperm,l_ri);
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
		printf("Step\t%d/%d\n",step_i,Nloop);
		addData1(XR,NX,folder_out,timestamp,"t_rad.dat");
		addData1(XZ,NX,folder_out,timestamp,"t_z.dat");
		addData1(XZ,NX,folder_out,timestamp,"t_tor.dat");
		

		

	//printf("ion:  0.\t %lf\t %lf\t %lf\n",XR[0],XZ[0],XT[0]);
	//printf("ion:  1.\t %lf\t %lf\t %lf\n",XR[1],XZ[1],XT[1]);
	//printf("ion:  0.\t %lf\t %lf\t %lf\n",VR[0],VZ[0],VT[0]);
	}	
	
	// Get CUDA timer 
	cudaEventElapsedTime(&time, start, stop);
	printf ("Time for the kernel: %f s\n", time/1000.0);

/*	// ION COORDS (device2HOST)
	cudaMemcpy(XR, xr, dimX, cudaMemcpyDeviceToHost);
	cudaMemcpy(XZ, xz, dimX, cudaMemcpyDeviceToHost);
	cudaMemcpy(XT, xt, dimX, cudaMemcpyDeviceToHost);

	// ION SPEEDS (device2HOST)
	cudaMemcpy(VR, vr, dimX, cudaMemcpyDeviceToHost);
	cudaMemcpy(VZ, vz, dimX, cudaMemcpyDeviceToHost);
	cudaMemcpy(VT, vt, dimX, cudaMemcpyDeviceToHost);*/

	//! MEMCOPY (device2HOST)
	cudaMemcpy(TMP, tmp, dimR, cudaMemcpyDeviceToHost);
	if(TMP[0]!=42.24){
		printf("\n+---	-------------------+\n | Fatal error in running. | \n | The CUDA did not run well. |\n+-----------------------+\n");
	}else{
		printf("\n	Memcopy OK.\n");
	}
	




	//! CUDA profiler STOP
	cudaProfilerStop();
	/*
	printf("ion:  0.\t %18.18le\t %18.18le\t %18.18le\n",VR[0],VZ[0],VT[0]);
	printf("ion:  0.\t %18.18le\t %18.18le\t %18.18le\n",XR[0],XZ[0],XT[0]);
	printf("----------------------------------------------------------\n");
*/

	printf("----------------------------------------------------------\n");
	printf("ion:  0.\t %lf\t %lf\t %lf\n",XR[0],XZ[0],XT[0]);
	printf("----------------------------------------------------------\n");
	for(int i=1; i<20; i++){
		printf("ion: %2d.\t %le\t %le\t %le\n",i,XR[i],XZ[i],XT[i]);
	}
	printf("----------------------------------------------------------\n");
/*
	//printf("ion:  0.\t %18.18le\t %18.18le\t %18.18le\n",XR[0],XZ[0],XT[0]);
	//printf("ion:  0.\t %18.18le\t %18.18le\t %18.18le\n",VR[0],VZ[0],VT[0]);
	
	*/
	//! Save data to files
	saveData1(XR,NX,folder_out,timestamp,"rad.dat");
	saveData1(XZ,NX,folder_out,timestamp,"z.dat");
	saveData1(XT,NX,folder_out,timestamp,"tor.dat");
	
	
	
	saveDataHT(concat("Shot ID: ",shotname),folder_out,timestamp);
	saveDataHT(concat("Run ID:  ",timestamp),folder_out,timestamp);
	saveDataHT("-----------------------------------",folder_out,timestamp);
	if(BANANA){
		saveDataHT("BANANA ORBITS",folder_out,timestamp);
	}else{		
		saveDataHT("ABP ION TRAJECTORIES",folder_out,timestamp);
		if(RADIONS){
			saveDataHT("(Real ionization position)",folder_out,timestamp);
			if($RENATE==110){
				saveDataHT("(TS + Renate 1.1.0)",folder_out,timestamp);
			}
		}else{
			saveDataHT("(R=const ionization)",folder_out,timestamp);
		}
	}
	saveDataHT("-----------------------------------",folder_out,timestamp);
	saveDataH("Beam energy","keV",energy,folder_out,timestamp);
	saveDataH("Atomic mass","AMU",mass,folder_out,timestamp);
	saveDataH("Beam diameter","mm",diameter,folder_out,timestamp);
	saveDataH2("Deflation (H/V)","°",$deflH,$deflV,folder_out,timestamp);
	if(!RADIONS&&!BANANA){	
		saveDataH("Ion. position (R)","m",R_midions,folder_out,timestamp);
	}
	
	saveDataH("Number of ions","",NX,folder_out,timestamp);
	saveDataHT("-----------------------------------",folder_out,timestamp);
	
	saveDataH("Detector position (R)","m",l_ri,folder_out,timestamp);
	
	saveDataHT("-----------------------------------",folder_out,timestamp);
	
	saveDataH("Timestep","s",dt,folder_out,timestamp);
	
	
	saveDataHT("-----------------------------------",folder_out,timestamp);
	
	saveDataH("Kernel runtime", "s", time/1000.0,folder_out,timestamp);
	saveDataHT("-----------------------------------",folder_out,timestamp);
	saveDataH("Number of blocks (threads)", "", max_blocks,folder_out,timestamp);
	saveDataH("Block size", "", BLOCK_SIZE,folder_out,timestamp);
	saveDataH("Length of a loop", "", Nstep,folder_out,timestamp);
	saveDataH("Number of loops", "", Nloop,folder_out,timestamp);
	
	

	printf("\nData folder: %s/%s\n\n",folder_out,timestamp);



	//! Free CUDA

	cudaFree(x_ptr);	cudaFree(xr);	cudaFree(xz);	cudaFree(xt);

	cudaFree(g_ptr);	cudaFree(rg);	cudaFree(zg);	
	
	cudaFree(br_ptr);	cudaFree(bz_ptr);	cudaFree(bt_ptr);	
	
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
	cudaFree(bt12);	cudaFree(bt13);	cudaFree(bt14);	cudaFree(bt15);	


	//! Free RAM
	free(RG);	free(ZG);	


	free(XR);	free(XZ);	free(XT);
	//	free(G_PTR);
	//	free(BR_PTR);	free(BZ_PTR);	free(BT_PTR);	
	
	free(BR0);	free(BR1);	free(BR2);	free(BR3);
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
	free(BT12);	free(BT13);	free(BT14);	free(BT15);	
	
	
	//! FREE TMP variables (RAM, cuda)
	free(TMP);	cudaFree(tmp);

	printf("Ready.\n\n");
}

char* concat(const char *s1, const char *s2){
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+1);//+1 for the zero-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}
