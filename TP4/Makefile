CC=gcc
CFLAGS=-Wall -O2 -std=c99 
LDLIBS=-pthread

all: build

build: main.c utils.c
	@$(CC) $(CFLAGS) $(LDLIBS) -o main $^
	
clean:
	rm -rf *.o

mrproper: clean
	rm -rf main
