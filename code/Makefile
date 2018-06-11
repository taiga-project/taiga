#nvcc -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32  main.cu -O3 -o taiga.exe

NVCC = nvcc
CFLAGS = -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32 -O3

SVNDEV = -D'SVN_REV="$(shell svnversion -n .)"' 

INPUT_FLAGS = $(SVNDEV)

all: taiga.exe taiga_renate.exe

taiga.exe: main.cu
	INPUT_FLAGS =  $(INPUT_FLAGS) -D'RENATE=0' 3DINPUTPROF=1 
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) -o taiga.exe  main.cu
	
	
taiga_renate.exe: main.cu
	INPUT_FLAGS = $(INPUT_FLAGS) -D'RENATE=110' 3DINPUTPROF=0
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) -o taiga_renate.exe  main.cu
	
clean: 
	rm taiga.exe
