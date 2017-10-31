#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "util.h"
#include "perf.h"

#define MIN 500
#define MAX 100000000
#define INNER 100
#define COEF 1.25
//1.25

long compute_ddot(long size, perf_t* start, perf_t* stop) {
  int i;
  double *dx = d_vector_create(size);
  double *dy = d_vector_create(size);

  d_vector_fill(dx, size, 2);
  d_vector_fill(dy, size, 5);

  perf(start);
  for (i = 0 ; i < INNER ; ++i) {
    cblas_ddot(size, dx, 1, dy, 1); 
  }
  perf(stop);

  d_vector_free(dx);
  d_vector_free(dy);  

  return (size + size - 1) * INNER;
}

int main(int argc, char** argv) {
  perf_t start, stop;
  double performance;
  long flop;
  long size = MIN;
  FILE* file = fopen("res.dot", "w");

  while (size < MAX) {  
    flop = compute_ddot(size, &start, &stop);

    perf_diff(&start, &stop);
    performance = perf_mflops(&stop, flop);
    printf("Size : %ld - Mflop/s : %lf \n", size, performance);
    fprintf(file, "%ld %lf\n", size, performance);
    size *= COEF;
  }
  fclose(file);
  return 0;
}

