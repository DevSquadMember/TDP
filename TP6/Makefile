MAKEFLAGS += -s
ifeq (run_seq,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

all: seq omp

seq:
	@gcc life_seq.c -o main_seq

omp:
	@g++ -fopenmp life_omp.c -o main_omp

run_seq: seq
	@./main_seq $(RUN_ARGS)

run_omp: omp
	@export OMP_NUM_THREADS=${NP}
	@./main_omp $(RUN_ARGS)

clean :
	@rm -rf *.o

mrproper: clean
	@rm -rf main_seq