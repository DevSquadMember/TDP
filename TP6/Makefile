MAKEFLAGS += -s
BUILD=build

ifeq (run_seq,$(firstword $(MAKECMDGOALS)))
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (run_omp,$(firstword $(MAKECMDGOALS)))
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (run_threads,$(firstword $(MAKECMDGOALS)))
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (run_mpi_sync,$(firstword $(MAKECMDGOALS)))
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (run_mpi_async,$(firstword $(MAKECMDGOALS)))
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (run_mpi_pers,$(firstword $(MAKECMDGOALS)))
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(RUN_ARGS):;@:)
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	C_OMP = clang-omp
else
	C_OMP = gcc
endif

all: dir seq omp threads mpi_sync

$(BUILD)/utils.o : utils.c utils.h
	@gcc -c utils.c -o $(BUILD)/utils.o

$(BUILD)/utils-mpi.o : utils-mpi.c utils-mpi.h
	@gcc -c utils-mpi.c -o $(BUILD)/utils-mpi.o

$(BUILD)/barrier.o : barrier.c barrier.h
	@gcc -c barrier.c -o $(BUILD)/barrier.o

dir: $(BUILD)

$(BUILD):
	@mkdir $(BUILD)

seq:
	@gcc life_seq.c -o $(BUILD)/main_seq

omp:
	@$(C_OMP) -fopenmp life_omp.c -o $(BUILD)/main_omp

threads: $(BUILD)/utils.o $(BUILD)/barrier.o
	@gcc -Wno-deprecated $(BUILD)/utils.o $(BUILD)/barrier.o life_threads.c -o $(BUILD)/main_threads

mpi_sync: $(BUILD)/utils.o $(BUILD)/barrier.o $(BUILD)/utils-mpi.o
	@mpicc -g -lm $(BUILD)/utils.o $(BUILD)/barrier.o $(BUILD)/utils-mpi.o life_mpi_sync.c -o $(BUILD)/main_mpi_sync

run_seq: seq
	@./$(BUILD)/main_seq $(RUN_ARGS)

run_omp: omp
	@env OMP_NUM_THREADS=${NP} ./$(BUILD)/main_omp $(RUN_ARGS)

run_threads: threads
	@./$(BUILD)/main_threads $(RUN_ARGS)

run_mpi_sync: mpi_sync
	@mpirun -np ${NP} $(BUILD)/main_mpi_sync $(RUN_ARGS)

clean :
	@rm -rf $(BUILD)/*.o

mrproper: clean
	@rm -rf $(BUILD)/*