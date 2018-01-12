MAKEFLAGS += -s

ifeq (run_seq,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (run_omp,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (run_threads,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

all: seq omp threads

seq:
	@gcc life_seq.c -o main_seq

omp:
	@gcc -fopenmp life_omp.c -o main_omp

threads:
	@gcc life_threads.c -o main_threads

run_seq: seq
	@./main_seq $(RUN_ARGS)

run_omp: omp
	@env OMP_NUM_THREADS=${NP} ./main_omp $(RUN_ARGS)

run_threads: threads
	@./main_threads $(RUN_ARGS)

clean :
	@rm -rf *.o

mrproper: clean
	@rm -rf main_seq
	@rm -rf main_omp
	@rm -rf main_threads