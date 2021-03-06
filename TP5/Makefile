MAKEFLAGS += -s
ifeq (generate,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

MAKEFLAGS += -s
ifeq (run,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

ifndef (NP)
	NP=3
endif

all: seq par

parser.o: parser.c parser.h
	@gcc -c parser.c

physics.o: physics.c physics.h
	@gcc -c physics.c

saver.o: saver.c saver.h
	@gcc -c saver.c

simulation.o : simulation.c simulation.h
	@gcc -c simulation.c

utils.o : utils.c utils.h
	@gcc -c utils.c

generator.o : generator.c
	@gcc -c generator.c

perf.o : perf.c perf.h
	@gcc -c perf.c

seq: simulation.o parser.o physics.o saver.o perf.o utils.o
	@gcc simulation.o parser.o physics.o saver.o perf.o utils.o main_seq.c -o main_seq

par: simulation.o parser.o physics.o saver.o perf.o utils.o
	@mpicc simulation.c parser.c physics.c saver.c perf.o utils.o main.c -o main

run_seq: all
	@./main_seq $(RUN_ARGS)

run: par
	@mpirun -np ${NP} main $(RUN_ARGS)

generate: generator.o parser.o physics.o simulation.o saver.o
	@gcc simulation.o parser.o physics.o saver.o generator.o -o generator
	@./generator $(RUN_ARGS)

compare: par simulation.o parser.o physics.o saver.o
	@gcc simulation.o parser.o physics.o saver.o compare.c -o comp
	@./comp

clean :
	@rm -rf *.o

mrproper: clean
	@rm -rf main
	@rm -rf main_par
	@rm -rf main_seq
	@rm -rf generator
	@rm -rf comp
	@rm -rf res.*