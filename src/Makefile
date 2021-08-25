NVCC = nvcc
COMMON_FLAGS = --ptxas-options=-v -w -Xptxas -dlcm=cg -maxrregcount=32 -lcurand -arch=compute_30 -I.
CFLAGS = $(COMMON_FLAGS) -O3
DEBUG_FLAGS = $(COMMON_FLAGS) -g -G
LIBS = -lcurand

VERSION = -D'TAIGA_VERSION="$(shell git branch | grep \* | cut -d ' ' -f2)"' -D'GIT_REV="$(shell git show -s --pretty=format:%h)"'

DEFAULT_FLAGS = $(VERSION) -D'RENATE=0' -D'READINPUTPROF=1' -D'FASTMODE=0'
RENATE_FLAGS = $(VERSION) -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=0'
RENATE_FAST_FLAGS = $(VERSION) -D'RENATE=110' -D'READINPUTPROF=0' -D'FASTMODE=1'

ifndef TAIGA_BUILD
	TAIGA_BUILD=build
endif
$(info Build directory: $(TAIGA_BUILD) )

all: taiga.exe taiga_debug.exe taiga_renate.exe taiga_renate_fast.exe t test_init field

no_test: taiga.exe taiga_debug.exe taiga_renate.exe taiga_renate_fast.exe

taiga.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/taiga.exe main.cu

taiga_debug.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(DEBUG_FLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/taiga_debug.exe main.cu

taiga_renate.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(RENATE_FLAGS) -o $(TAIGA_BUILD)/taiga_renate.exe main.cu

taiga_renate_fast.exe: main.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(RENATE_FAST_FLAGS) -o $(TAIGA_BUILD)/taiga_renate_fast.exe main.cu

t: $(TAIGA_BUILD) test

OBJ = $(TAIGA_BUILD)/obj

test: $(OBJ)/tests.o  $(OBJ)/test_bspline.o  $(OBJ)/test_solver.o $(OBJ)/test_basic_functions.o $(OBJ)/basic_functions.o
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) $^ -o $(TAIGA_BUILD)/test.exe

$(OBJ)/%.o: tests/%.c $(OBJ)
	$(NVCC)  $(CFLAGS) $(DEFAULT_FLAGS)  -I$tests -c $< -o $@

$(OBJ)/basic_functions.o: utils/basic_functions.c $(OBJ)
	$(NVCC)  $(CFLAGS) $(DEFAULT_FLAGS)  -I$utils -c $< -o $@

test_init: tests/test_taiga_init.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/test_init.exe tests/test_taiga_init.cu

test_framework: tests/taiga_test_example.c  | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/test_framework.exe tests/taiga_test_example.c

field: tests/test_field.cu | $(TAIGA_BUILD)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(TAIGA_BUILD)/test_field.exe tests/test_field.cu

$(TAIGA_BUILD):
	mkdir $@

$(OBJ):
	mkdir $@

$(CLEAN_O):
	rm $(TAIGA_BUILD)/$^

clean:
	rm $(TAIGA_BUILD)/*.exe
