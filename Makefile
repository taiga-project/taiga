NVCC = nvcc
GCC = gcc
COMMON_FLAGS = --ptxas-options=-v -w -Xptxas -dlcm=cg -maxrregcount=32 -lcurand -arch=compute_30 -Isrc -Iinclude
CFLAGS = $(COMMON_FLAGS) -O3
DEBUG_FLAGS = $(COMMON_FLAGS) -g -G
LIBS = -lcurand

VERSION = -D'TAIGA_VERSION="$(shell git branch | grep \* | cut -d ' ' -f2)"' -D'GIT_REV="$(shell git show -s --pretty=format:%h)"'

DEFAULT_FLAGS = $(VERSION) -D'RENATE=0' -D'FASTMODE=0'
RENATE_FLAGS = $(VERSION) -D'RENATE=1' -D'FASTMODE=0'
RENATE_FAST_FLAGS = $(VERSION) -D'RENATE=1' -D'FASTMODE=1'
TEST_FLAGS = $(VERSION) -fcommon

OBJ=build
BIN=bin

all: taiga.exe taiga_debug.exe taiga_renate.exe taiga_renate_fast.exe t test_init field

no_test: taiga.exe taiga_debug.exe taiga_renate.exe taiga_renate_fast.exe

taiga.exe: src/main.cu | $(BIN)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(BIN)/taiga.exe src/main.cu

taiga_debug.exe: src/main.cu | $(BIN)
	$(NVCC) $(DEBUG_FLAGS) $(DEFAULT_FLAGS) -o $(BIN)/taiga_debug.exe src/main.cu

taiga_renate.exe: src/main.cu | $(BIN)
	$(NVCC) $(CFLAGS) $(RENATE_FLAGS) -o $(BIN)/taiga_renate.exe src/main.cu

taiga_renate_fast.exe: src/main.cu | $(BIN)
	$(NVCC) $(CFLAGS) $(RENATE_FAST_FLAGS) -o $(BIN)/taiga_renate_fast.exe src/main.cu

t: test

test: $(OBJ)/tests.o  $(OBJ)/test_bspline.o  $(OBJ)/test_solver.o $(OBJ)/test_basic_functions.o $(OBJ)/basic_functions.o | $(BIN)
	$(GCC) $(TEST_FLAGS) -Isrc -Itests $^ -lm -o $(BIN)/test.exe

$(OBJ)/%.o: tests/%.c $(OBJ)
	$(GCC) $(TEST_FLAGS) -w -Isrc -Iinclude -Itests -c $< -lm -o $@

$(OBJ)/basic_functions.o: src/utils/basic_functions.c $(OBJ)
	$(GCC) $(TEST_FLAGS) -w -Isrc -Iinclude -Itests -I$utils -c $< -o $@

test_init: tests | $(BIN)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(BIN)/test_init.exe tests/test_taiga_init.cu

test_framework: tests  | $(BIN)
	$(GCC) $(TEST_FLAGS) -o $(BIN)/test_framework.exe tests/taiga_test_example.c

example_solvers: $(OBJ)/example_solvers.o $(OBJ)/test_solver.o | $(BIN)
	$(GCC) $(TEST_FLAGS) $^ -lm -o $(BIN)/example_solvers.exe

$(OBJ)/example_solvers.o: example/solvers/export_solver.c $(OBJ)
	$(GCC) $(TEST_FLAGS) -w -Isrc -Iinclude -Itests -c $< -lm -o $@

field: tests | $(BIN)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(BIN)/test_field.exe tests/test_field.cu

export_solver: tests | $(BIN)
	$(NVCC) $(CFLAGS) $(DEFAULT_FLAGS) -o $(BIN)/export_solver.exe tests/test_solver.c tests/export_solver.c

$(OBJ):
	mkdir $@

$(BIN):
	mkdir $@

$(CLEAN_O):
	rm $(OBJ)/$^

clean:
	rm bin/*
	rm build/*
