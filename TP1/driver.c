#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cblas.h"
#include "util.h"
#include "perf.h"

#define MIN 50
#define MAX 2000
#define INNER 100
#define COEF 1.25
#define ADD 0

void cblas_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  int n, m, k, val, b_pos, a_pos, c_pos;
  if (TransA == CblasNoTrans) {
    if (lda == 1) { 
      for (m = 0 ; m < M ; m++) {
        for (n = 0 ; n < N ; n++) {
          c_pos = n*ldc + m;
          val = C[c_pos];
          b_pos = n*ldb;
          for (k = 0 ; k < K ; k++) {
            val += A[k + m] * B[b_pos + k];
          } 
          C[c_pos] = alpha*val+beta;
        }
      }
    } else {
      for (m = 0 ; m < M ; m++) {
        for (n = 0 ; n < N ; n++) {
          c_pos = n*ldc + m;
          val = C[c_pos];
          b_pos = n*ldb;
          for (k = 0 ; k < K ; k++) {
            val += A[k*lda + m] * B[b_pos + k];
          } 
          C[c_pos] = alpha*val+beta;
        }
      }
    }
  } else if (TransA == CblasTrans) {
    if (alpha == 1 && beta == 0) {
      for (k = 0 ; k < K ; k++) {
        for (m = 0 ; m < M ; m++) {
          a_pos = k*lda + m;
          for (n = 0 ; n < N ; n++) {
            C[n*ldc + m] += A[a_pos] * B[n*ldb + k];
          } 
        }
      }
    } else if (alpha == 1) {
      for (k = 0 ; k < K ; k++) {
        for (m = 0 ; m < M ; m++) {
          a_pos = k*lda + m;
          for (n = 0 ; n < N ; n++) {
            C[n*ldc + m] += A[a_pos] * B[n*ldb + k] + beta;
          } 
        }
      }
    } else if (beta == 0) {
      for (k = 0 ; k < K ; k++) {
        for (m = 0 ; m < M ; m++) {
          a_pos = k*lda + m;
          for (n = 0 ; n < N ; n++) {
            C[n*ldc + m] += alpha*A[a_pos] * B[n*ldb + k];
          } 
        }
      }
    } else {
      for (k = 0 ; k < K ; k++) {
        for (m = 0 ; m < M ; m++) {
          a_pos = k*lda + m;
          for (n = 0 ; n < N ; n++) {
            C[n*ldc + m] += alpha*(A[a_pos] * B[n*ldb + k]) + beta;
          } 
        }
      }
    }
  } else {
    if (lda == 1) { 
      for (n = 0 ; n < N ; n++) {
        b_pos = n*ldb;
        a_pos = n*ldc;
        for (m = 0 ; m < M ; m++) {
          c_pos = a_pos + m;
          val = C[c_pos];
          for (k = 0 ; k < K ; k++) {
            val += alpha * A[k + m] * B[b_pos + k] + beta;
          } 
          C[c_pos] = val;
        }
      }
    } else {
      for (n = 0 ; n < N ; n++) {
        b_pos = n*ldb;
        a_pos = n*ldc;
        for (m = 0 ; m < M ; m++) {
          c_pos = a_pos + m;
          val = C[c_pos];
          for (k = 0 ; k < K ; k++) {
            val += alpha * A[k*lda + m] * B[b_pos + k] + beta;
          } 
          C[c_pos] = val;
        }
      }
    }
  }
}

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
  double *A = d_matrix_create(size, size);
  double *B = d_matrix_create(size, size);
  double *C = d_matrix_create(size, size);

  d_matrix_fill(A, size, size, 2);
  d_matrix_fill(B, size, size, 5);

  perf(start);
  cblas_dgemm_scalaire(CblasColMajor, trans, CblasNoTrans, size, size, size, 1, A, size, B, size, 0, C, size);
  perf(stop);

  d_matrix_free(A);
  d_matrix_free(B);
  d_matrix_free(C);

  if (trans == CblasNoTrans) {
    return (size * size * (5*size + 5));
  } else if (trans == CblasTrans) {
    return (size * size * (6*size + 1));
  } else {
    return (size * (2 + size * (7 * size + 1)));
  }
}

long compute_gemm_block(long size, perf_t* start, perf_t* stop) {
  double *A = d_matrix_create(size, size);
  double *B = d_matrix_create(size, size);
  double *C = d_matrix_create(size, size);

  d_matrix_fill(A, size, size, 2);
  d_matrix_fill(B, size, size, 5);
  d_matrix_fill(C, size, size, 0);

  perf(start);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, size, size, size, 1, A, size, B, size, 0, C, size);
  perf(stop);

  d_matrix_free(A);
  d_matrix_free(B);
  d_matrix_free(C);

  int nb = ceil(1.0*size/100);
  if (size <= 100) {
    return size*(size*(size*6 + 1));
  } else {
    return (((nb*5 + 3)*nb+1)*nb + size*(size*(size*6+1)));
  }
}

int main(int argc, char** argv) {
  perf_t start, stop;
  double performance, perf1, perf2, perf3;
  long flop;
  long size = MIN;
  FILE* file = fopen("res.dot", "w");

  while (size < MAX) { 

    /* DDOT */

    /*flop = compute_ddot(size, &start, &stop);
 *     perf_diff(&start, &stop);
 *         performance = perf_mflops(&stop, flop);
 *
 *             printf("Size : %ld - Mflops/s : %lf\n", size, performance);
 *                 fprintf(file, "%ld %lf\n", size, performance);
 *                 */
    /* GEMM */

    /*flop = compute_gemm(size, &start, &stop, CblasNoTrans);
 *     perf_diff(&start, &stop);
 *         perf1 = perf_mflops(&stop, flop);
 *
 *             flop = compute_gemm(size, &start, &stop, CblasTrans);
 *                 perf_diff(&start, &stop);
 *                     perf2 = perf_mflops(&stop, flop);
 *
 *                         flop = compute_gemm(size, &start, &stop, CblasConjTrans);
 *                             perf_diff(&start, &stop);
 *                                 perf3 = perf_mflops(&stop, flop);
 *
 *                                     printf("Size : %ld - Mflop/s : %lf - %lf - %lf\n", size, perf1, perf2, perf3);
 *                                         fprintf(file, "%ld %lf %lf %lf\n", size, perf1, perf2, perf3);
 *                                             */

    /* GEMM BLOCK */
    
    /*flop = compute_gemm_block(size, &start, &stop);
 *     perf_diff(&start, &stop);
 *         performance = perf_mflops(&stop, flop);
 *
 *             printf("Size : %ld - Mflop/s : %lf\n", size, performance);
 *                 fprintf(file, "%ld %lf\n", size, performance);*/

    /* GEMM VS GEMM BLOCK */

    flop = compute_gemm(size, &start, &stop, CblasNoTrans);
    perf_diff(&start, &stop);
    perf1 = perf_mflops(&stop, flop);

    flop = compute_gemm(size, &start, &stop, CblasTrans);
    perf_diff(&start, &stop);
    perf2 = perf_mflops(&stop, flop);

    flop = compute_gemm(size, &start, &stop, CblasConjTrans);
    perf_diff(&start, &stop);
    perf3 = perf_mflops(&stop, flop);

    flop = compute_gemm_block(size, &start, &stop);
    perf_diff(&start, &stop);
    performance = perf_mflops(&stop, flop);

    printf("Size : %ld - Mflop/s : %lf - %lf - %lf - %lf\n", size, perf1, perf2, perf3, performance);
    fprintf(file, "%ld %lf %lf %lf %lf\n", size, perf1, perf2, perf3, performance);
    
    size = size * COEF + ADD;
  }
  fclose(file);
  return 0;
}
