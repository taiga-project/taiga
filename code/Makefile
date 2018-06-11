#nvcc -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32  main.cu -O3 -o taiga.exe

NVCC = nvcc
CFLAGS = -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32 -O3

SVNDEV = -D'SVN_REV="$(shell svnversion -n .)"' 

INPUT_FLAGS = $(SVNDEV)
DEFAULT_FLAGS = -D'RENATE=0' -D'3DINPUTPROF=1'
RENATE_FLAGS = -D'RENATE=110' -D'3DINPUTPROF=0'

all: taiga.exe taiga_renate.exe

taiga.exe: main.cu
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) $(DEFAULT_FLAGS) -o taiga.exe  main.cu
	
	
taiga_renate.exe: main.cu
	
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) $(RENATE_FLAGS) -o taiga_renate.exe  main.cu
	
clean: 
	rm taiga.exe
