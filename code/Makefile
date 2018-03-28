#nvcc -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32  main.cu -O3 -o taiga.exe

NVCC = nvcc
CFLAGS = -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32 -O3
all: taiga.exe

taiga.exe: main.cu
	$(NVCC) $(CFLAGS) -o taiga.exe  main.cu
	
clean: rm taiga.exe
