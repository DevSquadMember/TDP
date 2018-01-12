MAKEFLAGS += -s
ifeq (run_seq,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

all: seq

seq:
	@gcc -fopenmp life_seq.c -o main_seq

run_seq: seq
	@./main_seq $(RUN_ARGS)
