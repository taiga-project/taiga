NVCC = nvcc
CFLAGS = --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32 -O3  -lcurand -arch=compute_30
TESTFLAGS = --ptxas-options=-v  -w  -Xptxas -dlcm=cg -maxrregcount=32 -g -G  -lcurand -arch=compute_30

VERSION = -D'TAIGA_VERSION="$(shell git branch | grep \* | cut -d ' ' -f2)"' -D'GIT_REV="$(shell git show -s --pretty=format:%h)"' 

INPUT_FLAGS = $(VERSION)
DEFAULT_FLAGS = -D'RENATE=0' -D'READINPUTPROF=1' -D'FASTMODE=0'
RENATE_FLAGS = -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=0'
RENATE_FAST_FLAGS = -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=1'

ifndef TAIGA_BUILD
	TAIGA_BUILD=build
endif
$(info Build directory: $(TAIGA_BUILD) )

all: taiga.exe taiga_debug.exe taiga_renate.exe taiga_renate_fast.exe

taiga.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/taiga.exe main.cu

taiga_debug.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(TESTFLAGS) $(INPUT_FLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/taiga_debug.exe main.cu

taiga_renate.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) $(RENATE_FLAGS) -o $(TAIGA_BUILD)/taiga_renate.exe main.cu

taiga_renate_fast.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(INPUT_FLAGS) $(RENATE_FAST_FLAGS) -o $(TAIGA_BUILD)/taiga_renate_fast.exe main.cu

t: test

test: tests/tests.cu | $(TAIGA_BUILD)
	$(NVCC) $(TESTFLAGS) $(INPUT_FLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/test.exe tests/tests.cu

$(TAIGA_BUILD):
	mkdir $@

clean:
	rm $(TAIGA_BUILD)/*.exe
