NVCC = nvcc
COMMON_FLAGS = --ptxas-options=-v -w -Xptxas -dlcm=cg -maxrregcount=32 -lcurand -arch=compute_30 -I.
CFLAGS = $(COMMON_FLAGS) -O3
TESTFLAGS = $(COMMON_FLAGS) -g -G

VERSION = -D'TAIGA_VERSION="$(shell git branch | grep \* | cut -d ' ' -f2)"' -D'GIT_REV="$(shell git show -s --pretty=format:%h)"'

DEFAULT_FLAGS = $(VERSION) -D'RENATE=0' -D'READINPUTPROF=1' -D'FASTMODE=0'
RENATE_FLAGS = $(VERSION) -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=0'
RENATE_FAST_FLAGS = $(VERSION) -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=1'

ifndef TAIGA_BUILD
	TAIGA_BUILD=build
endif
$(info Build directory: $(TAIGA_BUILD) )

all: taiga.exe taiga_debug.exe taiga_renate.exe taiga_renate_fast.exe

taiga.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/taiga.exe main.cu

taiga_debug.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(TESTFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/taiga_debug.exe main.cu

taiga_renate.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(RENATE_FLAGS) -o $(TAIGA_BUILD)/taiga_renate.exe main.cu

taiga_renate_fast.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(RENATE_FAST_FLAGS) -o $(TAIGA_BUILD)/taiga_renate_fast.exe main.cu

t: test bspline

test: tests/tests.cu | $(TAIGA_BUILD)
	$(NVCC) $(TESTFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/test.exe tests/tests.cu

test_framework: tests/taiga_test_example.c  | $(TAIGA_BUILD)
	$(NVCC) $(TESTFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/test_framework.exe tests/taiga_test_example.c

field: tests/test_field.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/test_field.exe tests/test_field.cu

bspline: tests/test_bspline.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/test_bspline.exe tests/test_bspline.cu

solver: tests/test_solver.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/test_solver.exe tests/test_solver.cu

$(TAIGA_BUILD):
	mkdir $@

clean:
	rm $(TAIGA_BUILD)/*.exe
