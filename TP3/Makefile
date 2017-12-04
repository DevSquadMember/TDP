MAKEFLAGS += -s

ifeq (run,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (generate,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

CFLAGS=-Wall -O2 -std=c99

all: build run

build:
	@mpicc $(CFLAGS) -lm -lmkl_intel_lp64 -lmkl_sequential -lmkl_core main.c parser.c perf.c utils.c -o main_par

run: all
	@export MKL_NUM_THREADS=1
	@mpirun -np ${NP} main_par $(RUN_ARGS)

generate: ./generator
	@gcc generator.c -o generator
	@./generator $(RUN_ARGS)

clean:
	@echo "Suppression des fichiers compilés..."
	@rm -f main_par
	@rm -f generator