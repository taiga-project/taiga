#nvcc -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32  main.cu -O3 -o taiga.exe

NVCC = nvcc
CFLAGS = -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32 -O3

SVNDEV = -D'SVN_REV="$(shell svnversion -n .)"' 

INPUT_FLAGS = $(SVNDEV)

all: taiga.exe

taiga.exe: main.cu
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) -D'RENATE=0' -o taiga.exe  main.cu
	
	
taiga_renate.exe: main.cu
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) -D'RENATE=110' -o taiga.exe  main.cu
	
clean: 
	rm taiga.exe
