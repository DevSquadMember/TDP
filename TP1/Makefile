CC=icc
CFLAGS=-Wall -O2 -std=c99 
LDLIBS=-pthread -mkl

all: driver

driver: driver.c util.c perf.c ddot.c dgemm_thread.c
	@$(CC) $(CFLAGS) $(LDLIBS) -o $@ $^
	
clean:
	rm -rf *.o

mrproper: clean
	rm -rf driver
