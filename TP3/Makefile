ifeq (run_par,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

CFLAGS=-Wall -O2 -std=c99
LDLIBS=-mkl

all: par

par:
	@mpicc $(CFLAGS) -lcblas main.c parser.c -o main_par

run_par: all
	@mpirun -np ${NP} main_par $(RUN_ARGS)

run_scat:
	@mpicc scatter.c parser.c -o scatter
	@mpirun -np 4 scatter

clean:
	@echo "Suppression des fichiers compilés..."
	@rm -f main_par
	@rm -f scatter