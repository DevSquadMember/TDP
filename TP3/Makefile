ifeq (run_par,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

all: par

par:
	@mpicc main.c parser.c -o main_par

run_par: all
	@mpirun -np ${NP} main_par $(RUN_ARGS)

run_scat:
	@mpicc scatter.c -o scatter
	@mpirun -np 4 scatter