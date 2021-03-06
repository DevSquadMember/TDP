MAKEFLAGS += -s
ifeq (run_par,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

ifeq (run_seq,$(firstword $(MAKECMDGOALS)))
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

ifndef (NP)
	NP=3
endif

all: par seq

par:
	@mpicc simulation.c parser.c physics.c saver.c main.c -o main_par

parser.o: parser.c parser.h
	@gcc -c parser.c

physics.o: physics.c physics.h
	@gcc -c physics.c

saver.o: saver.c saver.h
	@gcc -c saver.c

simulation.o : simulation.c simulation.h
	@gcc -c simulation.c

seq: simulation.o parser.o physics.o saver.o
	@gcc simulation.o parser.o physics.o saver.o main_seq.c -o main_seq

run_seq: all
	@./main_seq $(RUN_ARGS)

run_par: all
	@mpirun -np ${NP} main_par $(RUN_ARGS)

generate: ./generator
	@gcc generator.c -o generator
	@./generator $(RUN_ARGS)

compare: par simulation.o parser.o physics.o saver.o
	@gcc simulation.o parser.o physics.o saver.o compare.c -o comp
	@./comp

clean :
	@rm -rf *.o

mrproper: clean
	@rm -rf main_par
	@rm -rf main_seq
	@rm -rf generator
	@rm -rf comp
	@rm -rf res.*