#nvcc -arch sm_21 --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32  main.cu -O3 -o taiga.exe

NVCC = nvcc
CFLAGS = --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32 -O3  -lcurand -arch compute_30

VERSION = -D'TAIGA_VERSION="$(shell git branch | grep \* | cut -d ' ' -f2)"' -D'GIT_REV="$(shell git show -s --pretty=format:%h)"' 

INPUT_FLAGS = $(VERSION)
DEFAULT_FLAGS = -D'RENATE=0' -D'READINPUTPROF=1' -D'FASTMODE=0'
RENATE_FLAGS = -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=0'
RENATE_FAST_FLAGS = -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=1'

all: taiga.exe taiga_renate.exe taiga_renate_fast.exe

taiga.exe: main.cu
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) $(DEFAULT_FLAGS) -o taiga.exe  main.cu

taiga_renate.exe: main.cu
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) $(RENATE_FLAGS) -o taiga_renate.exe  main.cu

taiga_renate_fast.exe: main.cu
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) $(RENATE_FAST_FLAGS) -o taiga_renate_fast.exe  main.cu

clean:
	rm *.exe
