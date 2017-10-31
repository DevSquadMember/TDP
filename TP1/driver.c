#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "util.h"
#include "perf.h"

#define MIN 100
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

long compute_gemm(long size, perf_t* start, perf_t* stop, const enum CBLAS_TRANSPOSE trans) {
  int i;
  double *A = d_matrix_create(size, size);
  double *B = d_matrix_create(size, size);
  double *C = d_matrix_create(size, size);

  d_matrix_fill(A, size, size, 2);
  d_matrix_fill(B, size, size, 5);

  perf(start);
  cblas_dgemm(CblasColMajor, trans, CblasNoTrans, size, size, size, 1, A, size, B, size, 0, C, size);
  perf(stop);

  d_matrix_free(A);
  d_matrix_free(B);
  d_matrix_free(C);

  if (trans == CblasNoTrans) {
    return (size * size * (4*size + 3));
  } else if (trans == CblasTrans) {
    return (size * size * (6*size + 1));
  } else {
    return (size * (2 + size * (5 * size + 1)));
  }
}

int main(int argc, char** argv) {
  perf_t start, stop;
  double performance, perf1, perf2, perf3;
  long flop;
  long size = MIN;
  FILE* file = fopen("res.dot", "w");

  while (size < MAX) { 

    // DDOT

    flop = compute_ddot(size, &start, &stop);
    perf_diff(&start, &stop);
    performance = perf_mflops(&stop, flop);

    printf("Size : %ld - Mflops/s : %lf\n", size, performance);
    fprintf(file, "%ld %lf\n", size, performance);

    // GEMM

    /*flop = compute_gemm(size, &start, &stop, CblasNoTrans);
    perf_diff(&start, &stop);
    perf1 = perf_mflops(&stop, flop);

    flop = compute_gemm(size, &start, &stop, CblasTrans);
    perf_diff(&start, &stop);
    perf2 = perf_mflops(&stop, flop);

    flop = compute_gemm(size, &start, &stop, CblasConjTrans);
    perf_diff(&start, &stop);
    perf3 = perf_mflops(&stop, flop);

    printf("Size : %ld - Mflop/s : %lf - %lf - %lf\n", size, perf1, perf2, perf3);
    fprintf(file, "%ld %lf %lf %lf\n", size, perf1, perf2, perf3);*/
    size *= COEF;
  }
  fclose(file);
  return 0;
}

