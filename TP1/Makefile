CC=icc
CFLAGS=-Wall -O2 -std=c99 

all: driver

driver: driver.c util.c perf.c ddot.c dgemm.c
	@$(CC) $(CFLAGS) -o $@ $^
	
clean:
	rm -rf *.o

mrproper: clean
	rm -rf driver
