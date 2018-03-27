#nvcc -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32  main.cu -O3 -o taiga.exe

taiga.exe: main.cu
	nvcc -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32 -O3 -o 
