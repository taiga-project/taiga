// TAIGA default parameters

#define N_BLOCKS	 1		//! @param N_BLOCKS number of blocks (max 1M)
#define BLOCK_SIZE 	 192 		//! @param BLOCK_SIZE size of blocks (max 192 on Geforce GTS450) (max 768 on Geforce GTS650Ti)

#define $R_defl	2.3			//! radial position of deflection plates in meter -> TOROIDAL DEFLECTION
#define $deflH	 0				//! @param $deflH horizontal deflection in rad (up--down)  
#define $deflV	 0				//! @param $deflV vertical deflection in rad (left--right) -> TOROIDAL DEFLECTION

#define $default_energy   60				//! @param energy in keV
#define $default_mass  7.016004558			//! @param atomic mass in amu

#define $default_diameter 25//4/*e-20*/	  //! @param diameter in mm

#define dt	   1e-9			//! @param dt timestep in seconds

#define ERRORCHECK() cErrorCheck(__FILE__, __LINE__)
#define PI 3.141592653589792346
#define ELEMENTARY_CHARGE 1.60217656535e-19
#define AMU 1.66053892173e-27

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <time.h>
#include <math.h>
#include <string.h>
//#include <filesystem>

#include <cuda_profiler_api.h>
#include "cuda/nvToolsExt.h"

#include "main.cuh"
#include "dataio/fieldIn.c"

#if READINPUTPROF == 1
	#include "dataio/beamInFull.c"
#elif RENATE == 110
	#include "dataio/beamInRenate110.c"
#else
	#include "dataio/beamIn.c"
#endif

#include "dataio/beamOut.c"

#include "running/rk4.cu"
#include "running/ipol.cu"
#include "running/cyl2tor.cu"
#include "running/traj.cu"
#include "running/ctrl.cu"

#include "dataio/detectorIn.c"
#include "running/detector_postproc.cu"

int input_init_taiga(int argc, char *argv[], shot_prop *shot, beam_prop *beam){
	int max_blocks;

	if (argc > 1)	shot->name = argv[1];	
	if (argc > 2)	shot->runnumber = atoi(argv[2]);
	if (argc > 3)	beam->matter = argv[3];
	if (argc > 4)	beam->energy = atof(argv[4]);
	if (argc > 5)	beam->vertical_deflation = atof(argv[5]);
	if (argc > 6)	beam->diameter = atof(argv[6]);  

	if (argc > 8)	max_blocks = atoi(argv[8])/shot->block_size+1; 
		else	max_blocks=shot->block_number;

	if (argc > 9) shot->electric_field_module = atof(argv[9]);

	if (argc > 10){ 
		shot->step_host = atof(argv[10]); 
		shot->step_device = 1;
	}
	if (argc > 11)	shot->step_device = atof(argv[11]); 	
	if (argc > 12)	shot->debug = atof(argv[12]); 

	beam->mass = get_mass(beam->matter);
	return max_blocks;   
}

int main(int argc, char *argv[]){
	//! @param shotname name of shot folder input folder (8714,11344,11347)	
	
	shot_prop shot;
	beam_prop beam;
	int max_blocks = input_init_taiga(argc, argv, &shot, &beam);

	size_t dimD = 5 * sizeof(double);
	double *DETECTOR, *detector;
	DETECTOR = (double *)malloc(dimD);	cudaMalloc((void **) &detector,  dimD); 
	
	if (argc > 7)	fill_detector(DETECTOR, argv[7]);

	printf("shotname: %s\n",shot.name);  
	printf("detector: [ %lf %lf %lf %lf %lf]\n", DETECTOR[0],DETECTOR[1],DETECTOR[2],DETECTOR[3],DETECTOR[4]);

	int NX = shot.block_size * max_blocks;

	if (READINPUTPROF == 1){
		double *XR;
		NX = vectorReader0(&XR, "input/manual_profile/rad.dat");
		max_blocks = NX / shot.block_size+1;
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


	//! position and velocity array allocation
	size_t dimX = shot.block_size * max_blocks * sizeof(double);
	
	XR = (double*)malloc(dimX);
	XZ = (double*)malloc(dimX);
	XT = (double*)malloc(dimX);

	VR = (double*)malloc(dimX);
	VZ = (double*)malloc(dimX);
	VT = (double*)malloc(dimX);

	// phys. constants
	double eperm = ELEMENTARY_CHARGE/ AMU/ beam.mass;

	beamIn(XR, XZ, XT, VR, VZ, VT, beam.energy, eperm, NX, shot.name, beam.diameter, beam.toroidal_deflation, beam.vertical_deflation);

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
	size_t dimR = NR * sizeof(double);
	cudaMalloc((void **) &rg,  dimR); 
	
	//! Z-grid points
	int NZ = vectorReader(&ZG, "input/fieldSpl", shot.name, "z.spline");
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
	if (shot.electric_field_module)	shot.electric_field_module = electric_field_read_and_init(shot, &er_ptr,&ez_ptr,&et_ptr, dimRZ);
	
	// detector cell id
	size_t dimRint = NX * sizeof(int);
	int *DETCELLID, *detcellid;
	DETCELLID = (int *)malloc(dimRint);	cudaMalloc((void **) &detcellid,  dimRint); 	
	
	// temporary test data
	size_t dimService = 10 * sizeof(double);
	double *SERVICE_VAR, *service_var;
	SERVICE_VAR = (double *)malloc(dimService);	cudaMalloc((void **) &service_var,  dimService); 

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

	//! DETECTOR COORDS (HOST2device)
	cudaMemcpy(detector, DETECTOR, dimD, cudaMemcpyHostToDevice);
	
	// OUTPUT INIT
	addData1(XR,NX,folder_out,timestamp,"t_rad.dat");
	addData1(XZ,NX,folder_out,timestamp,"t_z.dat");
	addData1(XT,NX,folder_out,timestamp,"t_tor.dat");
	addData1(VR,NX,folder_out,timestamp,"t_vrad.dat");
	addData1(VZ,NX,folder_out,timestamp,"t_vz.dat");
	addData1(VT,NX,folder_out,timestamp,"t_vtor.dat");

	//! Set CUDA timer 
	cudaEvent_t start, stop;
	float runtime;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	if (shot.debug == 1)	debug_message_init(XR, XZ, XT, VR, VZ, VT);
	
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
				
		//ERRORCHECK();
		
		cudaEventRecord(start, 0);
		if (shot.electric_field_module){
			printf("electric_field_module ON\n");
			ctrl <<< max_blocks, shot.block_size >>> (NR,NZ,eperm,br_ptr,bz_ptr,bt_ptr,er_ptr,ez_ptr,et_ptr,g_ptr,x_ptr,v_ptr,detector,detcellid,shot.step_device,service_var);
		}else{
			ctrl <<< max_blocks, shot.block_size >>> (NR,NZ,eperm,br_ptr,bz_ptr,bt_ptr,g_ptr,x_ptr,v_ptr,detector,detcellid,shot.step_device,service_var);
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
		
		// DETCELLID (device2HOST)
		cudaMemcpy(DETCELLID, detcellid, dimRint, cudaMemcpyDeviceToHost);
		//ERRORCHECK();
		
		// Save data to files
		printf("Step\t%d/%d\n",step_i,shot.step_host);
		addData1(XR,NX,folder_out,timestamp,"t_rad.dat");
		addData1(XZ,NX,folder_out,timestamp,"t_z.dat");
		addData1(XT,NX,folder_out,timestamp,"t_tor.dat");
		addData1(VR,NX,folder_out,timestamp,"t_vrad.dat");
		addData1(VZ,NX,folder_out,timestamp,"t_vz.dat");
		addData1(VT,NX,folder_out,timestamp,"t_vtor.dat");
		
		if (shot.debug == 1)	debug_message_run(XR, XZ, XT, VR, VZ, VT);

		/*if (shot.step_host > 1){
			for (int i = 1; (i < NX && XR[i] == detector); i++){;
				if (i == NX-1) shot.step_host = step_i;
			}
		}*/
	}

	// Get CUDA timer 
	cudaEventElapsedTime(&runtime, start, stop);
	printf ("Time for the kernel: %f s\n", runtime/1000.0);

	//! MEMCOPY (device2HOST)
	cudaMemcpy(SERVICE_VAR, service_var, dimService, cudaMemcpyDeviceToHost);
	if(SERVICE_VAR[0]!=42.24){
		printf("\n +--------------------------+\n | Fatal error in running.    | \n | The CUDA did not run well. |\n +---------------------------+\n");
	}else{
		printf("\n	Memcopy OK.\n");
	}

	/*if (shot.debug == 1){
		for (int i=0;i<10;i++) {
			printf("SERVICE_VAR%d\t%lf\n",i,SERVICE_VAR[i]);
		}
	}*/

	detector_module(x_ptr, detector, detcellid, "test"); //detector_name = "test"

	//! CUDA profiler STOP
	cudaProfilerStop();

	//! Save data to files
	saveData1(XR,NX,folder_out,timestamp,"rad.dat");
	saveData1(XZ,NX,folder_out,timestamp,"z.dat");
	saveData1(XT,NX,folder_out,timestamp,"tor.dat");
	saveData1(VR,NX,folder_out,timestamp,"vrad.dat");
	saveData1(VZ,NX,folder_out,timestamp,"vz.dat");
	saveData1(VT,NX,folder_out,timestamp,"vtor.dat");
	saveData1(DETCELLID,NX,folder_out,timestamp,"detcellid.dat");
	
	saveDataHT(concat("Shot ID: ",shot.name),folder_out,timestamp);
	saveDataHT(concat("Run ID:  ",timestamp),folder_out,timestamp);
	saveDataHT("-----------------------------------",folder_out,timestamp);
	saveDataHT(concat("version: r ",SVN_REV),folder_out,timestamp);
	saveDataHT("-----------------------------------",folder_out,timestamp);
		
	saveDataHT("ABP ION TRAJECTORIES",folder_out,timestamp);

	saveDataHT("(Real ionization position)",folder_out,timestamp); 
	if(READINPUTPROF==1){
		saveDataHT("(3D input)",folder_out,timestamp);			
	}else if(RENATE==110){
		saveDataHT("(TS + Renate 1.1.0)",folder_out,timestamp);
	}

	saveDataHT("-----------------------------------",folder_out,timestamp);

	if(!READINPUTPROF){
		saveDataH("Beam energy","keV",beam.energy,folder_out,timestamp);
		saveDataH("Atomic mass","AMU",beam.mass,folder_out,timestamp);
		saveDataH("Beam diameter","mm",beam.diameter,folder_out,timestamp);
		saveDataH2("Deflation (toroidal/vertical)","°",beam.toroidal_deflation,beam.vertical_deflation,folder_out,timestamp);
	}
	
	
	saveDataH("Number of ions","",NX,folder_out,timestamp);
	saveDataHT("-----------------------------------",folder_out,timestamp);
	
	 // DETECTOR
	saveDataH("Detector position (R)","m",DETECTOR[0],folder_out,timestamp);
	saveDataH("Detector position (Z)","m",DETECTOR[1],folder_out,timestamp);
	saveDataH("Detector position (T)","m",DETECTOR[2],folder_out,timestamp);
	saveDataH("Detector angle (Z/R)","°",atan(DETECTOR[3])/PI*180.0,folder_out,timestamp);
	saveDataH("Detector angle (T/R)","°",atan(DETECTOR[4])/PI*180.0,folder_out,timestamp);
	
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

	//! Free RAM
	free(RG);	free(ZG);	
	free(XR);	free(XZ);	free(XT);
	
	//! FREE SERVICE_VAR variables (RAM, cuda)
	free(SERVICE_VAR);	cudaFree(service_var);

	printf("Ready.\n\n");
}


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
		mass = 22.98976928;
	}else if (strcmp(s,"K")==0){
		mass = 39.9639984821;
	}else if (strcmp(s,"H2")==0){
		mass = 2.013553212724;
	}else if (strcmp(s,"Li7")==0){
		mass = 7.016004558;
	}else if (strcmp(s,"Na23")==0){
		mass = 22.98976928;
	}else if (strcmp(s,"K40")==0){
		mass = 39.9639984821;
	}else{
		try{
			mass = atof(s);
		}catch (...){
			mass = (double)$default_mass;
		}
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
	
	if (shot.debug == 1){
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

void fill_detector(double *DETECTOR, char* values){

	char *el; 
	el = strtok(values,",");	DETECTOR[0] = strtod (el, NULL);
	el = strtok(NULL,",");	DETECTOR[1] = strtod (el, NULL);
	el = strtok(NULL,",");	DETECTOR[2] = tan(strtod (el, NULL) * PI/180.0);
	el = strtok(NULL,",");	DETECTOR[3] = tan(strtod (el, NULL) * PI/180.0);

}


void process_detector(int *detcellid, double **x_ptr){
	
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

//

char* concat(const char *s1, const char *s2){
	char *result = (char*)malloc(strlen(s1)+strlen(s2)+1);
	strcpy(result, s1);
	strcat(result, s2);
	return result;
}

// DEBUG

void debug_message_init(double* XR, double* XZ, double* XT, double* VR, double* VZ, double* VT){
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

void debug_message_run(double* XR, double* XZ, double* XT, double* VR, double* VZ, double* VT){
			printf("Xion:  0.\t %lf\t %lf\t %lf\n",XR[0],XZ[0],XT[0]);
			printf("Xion:  1.\t %lf\t %lf\t %lf\n",XR[1],XZ[1],XT[1]);
			printf("Vion:  0.\t %lf\t %lf\t %lf\n",VR[0],VZ[0],VT[0]);
}
