ifeq (run,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (check,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

CC=gcc
CFLAGS=-Wall -O2 -std=c99 
LDLIBS=-pthread

all: build

build: lib_utils.c lib_matrix.c timer.c perf.c utils.c
	@mpicc $(CFLAGS) $(LDLIBS) -o main $^ main.c
	@mpicc $(CFLAGS) $(LDLIBS) -o check $^ check.c

run: all
	@mpirun -np ${NP} main $(RUN_ARGS)

check: all
	@mpirun -np ${NP} check $(RUN_ARGS)
	
clean:
	rm -rf *.o

mrproper: clean
	rm -rf main
